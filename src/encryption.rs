use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, PrimeField};
use ark_serialize::*;
use ark_std::UniformRand;
use merlin::Transcript;
use rand::thread_rng;
use retry::{delay::NoDelay, retry};

use crate::utils::{add_to_transcript, hash_to_bytes, xor};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Default)]
pub struct DLogProof<F: PrimeField> {
    pub c: F,       //challenge
    pub z_alpha: F, //opening for alpha
    pub z_beta: F,  //opening for beta
    pub z_s: F,     //opening for s
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Ciphertext<E: Pairing> {
    pub ct1: [u8; 32],
    pub ct2: E::G2,
    pub ct3: E::G2,
    pub ct4: E::G2,
    pub gs: E::G1,
    pub x: E::ScalarField,
    pub pi: DLogProof<E::ScalarField>,
}

impl<E: Pairing> Ciphertext<E> {
    /// panicks if ciphertext does not verify
    pub fn verify(&self, htau: E::G2, pk: E::G2) {
        let g = E::G1::generator();
        let h = E::G2::generator();

        // k2.ct2^c = h^{(tau-x)*z_alpha}, k3.ct3^c = h^{z_alpha} * pk^{z_beta}, k4.ct4^c = h^{z_beta}, and k_s.gs^c = g^{z_s}
        let minus_c = -self.pi.c;
        let recovered_k2 = (htau - (h * self.x)) * self.pi.z_alpha + (self.ct2 * minus_c);
        let recovered_k3 = h * self.pi.z_alpha + pk * self.pi.z_beta + (self.ct3 * minus_c);
        let recovered_k4 = h * self.pi.z_beta + (self.ct4 * minus_c);
        let recovered_k_s = g * self.pi.z_s + (self.gs * minus_c);

        let mut ts: Transcript = Transcript::new(&[0u8]);
        add_to_transcript(&mut ts, b"ct1", self.ct1);
        add_to_transcript(&mut ts, b"ct2", self.ct2);
        add_to_transcript(&mut ts, b"ct3", self.ct3);
        add_to_transcript(&mut ts, b"ct4", self.ct4);
        add_to_transcript(&mut ts, b"gs", self.gs);
        add_to_transcript(&mut ts, b"x", self.x);

        add_to_transcript(&mut ts, b"k2", recovered_k2);
        add_to_transcript(&mut ts, b"k3", recovered_k3);
        add_to_transcript(&mut ts, b"k4", recovered_k4);
        add_to_transcript(&mut ts, b"k_s", recovered_k_s);

        // Fiat-Shamir to get challenge
        let mut c_bytes = [0u8; 31];
        ts.challenge_bytes(&[8u8], &mut c_bytes);
        let c = E::ScalarField::from_random_bytes(&c_bytes).unwrap();

        // assert that the recomputed challenge matches
        assert_eq!(self.pi.c, c);
    }
}

pub fn encrypt<E: Pairing>(
    msg: [u8; 32],
    x: E::ScalarField,
    hid: E::G1,
    htau: E::G2,
    pk: E::G2,
) -> Ciphertext<E> {
    let rng = &mut thread_rng();

    let g = E::G1::generator();
    let h = E::G2::generator();

    // hash element S to curve to get tg
    // retry if bytes cannot be converted to a field element
    let result = retry(NoDelay, || {
        let s = E::ScalarField::rand(rng);
        let gs = g * s;
        let hgs = hash_to_bytes(gs);
        let tg_option = E::ScalarField::from_random_bytes(&hgs);

        match tg_option {
            Some(tg) => Ok((s, gs, tg)),
            None => {
                #[cfg(debug_assertions)]
                {
                    dbg!("Failed to hash to field element, retrying...");
                }
                Err(())
            }
        }
    });

    let (s, gs, tg) = result.unwrap();

    // compute mask
    let alpha = E::ScalarField::rand(rng);
    let beta = E::ScalarField::rand(rng);
    let mask = E::pairing(hid - (g * tg), h) * alpha; //e(H(id)/g^tg, h)^alpha
    let hmask = hash_to_bytes(mask);

    // xor msg and hmask
    let ct1: [u8; 32] = xor(&msg, &hmask).as_slice().try_into().unwrap();
    let ct2 = (htau - (h * x)) * alpha; //h^{(tau-x)*alpha}
    let ct3 = h * alpha + pk * beta; //h^alpha * pk^beta
    let ct4 = h * beta; //h^beta

    // prove knowledge of alpha, beta, and s such that ct2  = h^{(tau-x)*alpha}, ct3 = h^alpha * pk^beta, ct4 = h^beta, and gs = g^s
    // prover sends k2 = h^{(tau-x)*r_alpha}, k3 = h^{r_alpha} * pk^{r_beta}, k4 = h^{r_beta}, and k_s = g^{r_s}
    // verifier sends a random challenge c
    // prover sends z_alpha = r_alpha + c*alpha, z_beta = r_beta + c*beta, and z_s = r_s + c*s
    // verifier checks that k2.ct2^c = h^{(tau-x)*z_alpha}, k3.ct3^c = h^{z_alpha} * pk^{z_beta}, k4.ct4^c = h^{z_beta}, and k_s.gs^c = g^{z_s}

    let r_alpha = E::ScalarField::rand(rng);
    let r_beta = E::ScalarField::rand(rng);
    let r_s = E::ScalarField::rand(rng);

    let k2 = (htau - (h * x)) * r_alpha;
    let k3 = h * r_alpha + pk * r_beta;
    let k4 = h * r_beta;
    let k_s = g * r_s;

    let mut ts: Transcript = Transcript::new(&[0u8]);
    add_to_transcript(&mut ts, b"ct1", ct1);
    add_to_transcript(&mut ts, b"ct2", ct2);
    add_to_transcript(&mut ts, b"ct3", ct3);
    add_to_transcript(&mut ts, b"ct4", ct4);
    add_to_transcript(&mut ts, b"gs", gs);
    add_to_transcript(&mut ts, b"x", x);

    add_to_transcript(&mut ts, b"k2", k2);
    add_to_transcript(&mut ts, b"k3", k3);
    add_to_transcript(&mut ts, b"k4", k4);
    add_to_transcript(&mut ts, b"k_s", k_s);

    // Fiat-Shamir to get challenge
    let mut c_bytes = [0u8; 31];
    ts.challenge_bytes(&[8u8], &mut c_bytes);
    let c = E::ScalarField::from_random_bytes(&c_bytes).unwrap();

    let z_alpha = r_alpha + c * alpha;
    let z_beta = r_beta + c * beta;
    let z_s = r_s + c * s;

    let pi = DLogProof {
        c,
        z_alpha,
        z_beta,
        z_s,
    };

    Ciphertext {
        ct1,
        ct2,
        ct3,
        ct4,
        gs,
        x,
        pi,
    }
}

#[cfg(test)]
mod tests {
    use crate::dealer::Dealer;

    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::bls12::Bls12;
    use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
    use rand::thread_rng;

    type E = Bls12_381;
    type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
    type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
    type G2 = <Bls12<ark_bls12_381::Config> as Pairing>::G2;

    #[test]
    fn test_encryption() {
        let mut rng = thread_rng();

        let batch_size = 1 << 5;
        let n = 1 << 4;
        let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();

        let mut dealer = Dealer::<E>::new(batch_size, n, n / 2 - 1);
        let (crs, _) = dealer.setup(&mut rng);
        let pk = dealer.get_pk();

        let msg = [1u8; 32];
        let x = tx_domain.group_gen;

        let hid = G1::rand(&mut rng);

        let ct = encrypt::<Bls12_381>(msg, x, hid, crs.htau, pk);

        let mut ct_bytes = Vec::new();
        ct.serialize_compressed(&mut ct_bytes).unwrap();
        println!("Compressed ciphertext: {} bytes", ct_bytes.len());

        let mut ct_bytes = Vec::new();
        ct.serialize_uncompressed(&mut ct_bytes).unwrap();
        println!("Uncompressed ciphertext: {} bytes", ct_bytes.len());

        let mut g1_bytes = Vec::new();
        let mut g2_bytes = Vec::new();
        let mut fr_bytes = Vec::new();

        let g = G1::generator();
        let h = G2::generator();
        let x = tx_domain.group_gen;

        g.serialize_compressed(&mut g1_bytes).unwrap();
        h.serialize_compressed(&mut g2_bytes).unwrap();
        x.serialize_compressed(&mut fr_bytes).unwrap();

        println!("G1 len: {} bytes", g1_bytes.len());
        println!("G2 len: {} bytes", g2_bytes.len());
        println!("Fr len: {} bytes", fr_bytes.len());

        ct.verify(crs.htau, pk);
    }
}
