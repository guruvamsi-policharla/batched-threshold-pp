use ark_ec::{pairing::Pairing, CurveGroup, Group, VariableBaseMSM};
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use ark_serialize::CanonicalSerialize;
use ark_std::{One, Zero};
use merlin::Transcript;
use std::ops::Div;

use crate::dealer::CRS;

pub fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(!v.is_empty());
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}

pub fn hash_to_bytes<T: CanonicalSerialize>(inp: T) -> [u8; 32] {
    let mut bytes = Vec::new();
    inp.serialize_uncompressed(&mut bytes).unwrap();
    let hash = blake3::hash(bytes.as_slice());
    let hash_bytes = hash.as_bytes();
    *hash_bytes
}

pub fn xor(a: &[u8], b: &[u8]) -> Vec<u8> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| x ^ y).collect()
}

pub fn add_to_transcript<T: CanonicalSerialize>(
    ts: &mut Transcript,
    label: &'static [u8],
    data: T,
) {
    let mut data_bytes = Vec::new();
    data.serialize_uncompressed(&mut data_bytes).unwrap();
    ts.append_message(label, &data_bytes);
}

/// compute KZG opening proof
pub fn compute_opening_proof<E: Pairing>(
    crs: &CRS<E>,
    polynomial: &DensePolynomial<E::ScalarField>,
    point: &E::ScalarField,
) -> E::G1 {
    let eval = polynomial.evaluate(point);
    let eval_as_poly = DensePolynomial::from_coefficients_vec(vec![eval]);
    let numerator = polynomial - &eval_as_poly;
    let divisor = DensePolynomial::from_coefficients_vec(vec![
        E::ScalarField::zero() - point,
        E::ScalarField::one(),
    ]);
    let witness_polynomial = numerator.div(&divisor);

    pedersen_commit::<E::G1>(&crs.powers_of_g, &witness_polynomial.coeffs())
}

pub fn pedersen_commit<G: Group + CurveGroup>(bases: &[G], scalars: &[G::ScalarField]) -> G {
    debug_assert_eq!(bases.len(), scalars.len());

    let plain_scalars = convert_to_bigints(&scalars);
    let affine_bases = bases.iter().map(|g| g.into_affine()).collect::<Vec<_>>();
    <G as VariableBaseMSM>::msm_bigint(&affine_bases, &plain_scalars)
}

fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    let coeffs = ark_std::cfg_iter!(p)
        .map(|s| s.into_bigint())
        .collect::<Vec<_>>();
    coeffs
}

/// Computes all the openings of a KZG commitment in O(n log n) time
/// See https://github.com/khovratovich/Kate/blob/master/Kate_amortized.pdf
/// eprint version has a bug and hasn't been updated
pub fn open_all_values<E: Pairing>(
    y: &Vec<E::G1>,
    f: &Vec<E::ScalarField>,
    domain: &Radix2EvaluationDomain<E::ScalarField>,
) -> Vec<E::G1> {
    let top_domain = Radix2EvaluationDomain::<E::ScalarField>::new(2 * domain.size()).unwrap();

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    // f = {f0 ,f1, ..., fd}
    // v = {(d 0s), f1, ..., fd}
    let mut v = vec![E::ScalarField::zero(); domain.size() + 1];
    v.append(&mut f[1..f.len()].to_vec());

    debug_assert_eq!(v.len(), 2 * domain.size());
    let v = top_domain.fft(&v);

    // h = y \odot v
    let mut h = vec![E::G1::zero(); 2 * domain.size()];
    for i in 0..2 * domain.size() {
        h[i] = y[i] * (v[i]);
    }

    // inverse fft on h
    let mut h = top_domain.ifft(&h);

    h.truncate(domain.size());

    // fft on h to get KZG proofs
    let pi = domain.fft(&h);

    pi
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::{bls12::Bls12, pairing::Pairing, Group};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{UniformRand, Zero};

    use crate::dealer::Dealer;

    use super::*;
    type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
    type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
    type G2 = <Bls12<ark_bls12_381::Config> as Pairing>::G2;
    type E = Bls12_381;

    #[test]
    fn open_all_test() {
        let mut rng = ark_std::test_rng();

        let domain_size = 1 << 5;
        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();

        let mut dealer = Dealer::<E>::new(domain_size, 1 << 5, domain_size / 2 - 1);
        let (crs, _) = dealer.setup(&mut rng);

        let mut f = vec![Fr::zero(); domain_size];
        for i in 0..domain_size {
            f[i] = Fr::rand(&mut rng);
        }

        let com = pedersen_commit::<G1>(&crs.powers_of_g, &f);
        let pi = open_all_values::<E>(&crs.y, &f, &domain);

        // verify the kzg proof
        let g = G1::generator();
        let h = G2::generator();

        let fpoly = DensePolynomial::from_coefficients_vec(f.clone());
        for i in 0..domain_size {
            let lhs = E::pairing(com - (g * fpoly.evaluate(&domain.element(i))), h);
            let rhs = E::pairing(pi[i], crs.htau - (h * domain.element(i)));
            assert_eq!(lhs, rhs);
        }
    }
}
