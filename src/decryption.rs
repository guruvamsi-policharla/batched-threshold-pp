use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::Zero;

use crate::{
    dealer::CRS,
    encryption::Ciphertext,
    utils::{hash_to_bytes, open_all_values, pedersen_commit, xor},
};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretKey<E: Pairing> {
    sk_share: E::ScalarField,
}

impl<E: Pairing> SecretKey<E> {
    pub fn new(sk_share: E::ScalarField) -> Self {
        SecretKey { sk_share }
    }

    /// each party in the committee computes a partial decryption
    pub fn partial_decrypt(
        &self,
        ct: &Vec<Ciphertext<E>>,
        hid: E::G1,
        pk: E::G2,
        crs: &CRS<E>,
    ) -> E::G1 {
        let batch_size = crs.powers_of_g.len();
        for i in 0..batch_size {
            ct[i].verify(crs.htau, pk);
        }

        let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(batch_size).unwrap();

        let mut fevals = vec![E::ScalarField::zero(); batch_size];
        for i in 0..batch_size {
            let tg_bytes = hash_to_bytes(ct[i].gs);
            fevals[i] = E::ScalarField::from_random_bytes(&tg_bytes).unwrap();
        }
        let fcoeffs = tx_domain.ifft(&fevals);
        let com = pedersen_commit::<E::G1>(&crs.powers_of_g, &fcoeffs);
        let delta = hid - com;

        let pd = delta * self.sk_share;

        pd
    }
}

/// decrypts all the ciphertexts in a batch
pub fn decrypt_all<E: Pairing>(
    partial_decryptions: &Vec<E::G1>,
    ct: &Vec<Ciphertext<E>>,
    hid: E::G1,
    crs: &CRS<E>,
) -> Vec<[u8; 32]> {
    let batch_size = ct.len();

    let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(batch_size).unwrap();

    // compute fevals by hashing gs of the ciphertexts to get fevals
    let mut fevals = vec![E::ScalarField::zero(); batch_size];
    for i in 0..batch_size {
        let tg_bytes = hash_to_bytes(ct[i].gs);
        fevals[i] = E::ScalarField::from_random_bytes(&tg_bytes).unwrap();
    }

    let fcoeffs = tx_domain.ifft(&fevals);

    let com = pedersen_commit::<E::G1>(&crs.powers_of_g, &fcoeffs);
    let delta = hid - com;

    // compute msm with lagrange_coeffs_0 and partial_decryptions
    let sigma = pedersen_commit::<E::G1>(&partial_decryptions, &crs.lagrange_coeffs_0);

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    let pi = open_all_values::<E>(&crs.y, &fcoeffs, &tx_domain);

    // now decrypt each of the ciphertexts as m = ct(1) xor H(e(delta,ct2).e(pi,ct3)e(-sigma,ct4))
    let mut m = vec![[0u8; 32]; batch_size];
    for i in 0..batch_size {
        let mask = E::pairing(pi[i], ct[i].ct2)
            + E::pairing(delta, ct[i].ct3)
            + E::pairing(-sigma, ct[i].ct4);
        let hmask = hash_to_bytes(mask);
        m[i] = xor(&ct[i].ct1, &hmask).as_slice().try_into().unwrap();
    }

    m
}
