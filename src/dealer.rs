use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, Group};
use ark_ff::PrimeField;
use ark_poly::{domain::EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{rand::RngCore, One, UniformRand, Zero};
use rand::thread_rng;
use std::{iter, vec};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct CRS<E: Pairing> {
    pub powers_of_g: Vec<E::G1>,
    pub htau: E::G2,

    pub y: Vec<E::G1>,

    pub lagrange_coeffs_0: Vec<E::ScalarField>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Dealer<E: Pairing> {
    batch_size: usize,
    n: usize,
    t: usize, // t+1 parties need to agree to decrypt
    sk: E::ScalarField,
}

impl<E> Dealer<E>
where
    E: Pairing,
{
    pub fn new(batch_size: usize, n: usize, t: usize) -> Self {
        let rng = &mut thread_rng();
        Self {
            batch_size,
            n,
            t,
            sk: E::ScalarField::rand(rng),
        }
    }

    pub fn get_pk(&self) -> E::G2 {
        E::G2::generator() * self.sk
    }

    pub fn setup<R: RngCore>(&mut self, rng: &mut R) -> (CRS<E>, Vec<E::ScalarField>) {
        // Sample tau and compute its powers ==========================================================
        let tau = E::ScalarField::rand(rng);
        let powers_of_tau: Vec<E::ScalarField> =
            iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(self.batch_size)
                .collect();

        // Generators
        let g = E::G1::generator();
        let h = E::G2::generator();

        // Compute powers of g
        let window_size = FixedBase::get_mul_window_size(self.batch_size);
        let scalar_size = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let g_table = FixedBase::get_window_table(scalar_size, window_size, g);
        let powers_of_g =
            FixedBase::msm::<E::G1>(scalar_size, window_size, &g_table, &powers_of_tau);

        // Compute the Toeplitz matrix preprocessing ==================================================
        let mut top_tau = powers_of_tau.clone();
        top_tau.truncate(self.batch_size);
        top_tau.reverse();
        top_tau.resize(2 * self.batch_size, E::ScalarField::zero());

        let top_domain =
            Radix2EvaluationDomain::<E::ScalarField>::new(2 * self.batch_size).unwrap();
        let top_tau = top_domain.fft(&top_tau);

        // Compute powers of top_tau
        let window_size = FixedBase::get_mul_window_size(2 * self.batch_size);
        let scalar_size = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let top_tau_table = FixedBase::get_window_table(scalar_size, window_size, g);
        let y = FixedBase::msm::<E::G1>(scalar_size, window_size, &top_tau_table, &top_tau);

        let mut sk_shares = vec![E::ScalarField::zero(); self.n];
        sk_shares[0] = self.sk;
        for i in 1..self.t {
            sk_shares[i] = E::ScalarField::rand(rng);
        }

        let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        share_domain.fft_in_place(&mut sk_shares);

        let lagrange_coeffs_0 =
            share_domain.evaluate_all_lagrange_coefficients(E::ScalarField::zero());

        let crs = CRS::<E> {
            powers_of_g,
            htau: h * tau,
            y,
            lagrange_coeffs_0,
        };

        (crs, sk_shares)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    type E = Bls12_381;
    type Fr = <E as Pairing>::ScalarField;
    type G2 = <E as Pairing>::G2;

    #[test]
    fn test_dealer() {
        let mut rng = ark_std::test_rng();
        let batch_size = 1 << 5;
        let n = 1 << 4;
        let t = n / 2 - 1;

        let mut dealer = Dealer::<E>::new(batch_size, n, t);
        let (crs, sk_shares) = dealer.setup(&mut rng);

        let share_domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
        let should_be_sk = share_domain.ifft(&sk_shares)[0];

        let pk = dealer.get_pk();
        let g_sks = sk_shares
            .iter()
            .map(|ski| G2::generator() * ski)
            .collect::<Vec<_>>();

        let should_be_pk = share_domain.ifft(&g_sks)[0];
        assert_eq!(pk, should_be_pk);

        assert_eq!(dealer.sk, should_be_sk);
        assert_eq!(crs.powers_of_g.len(), batch_size);
        assert_eq!(sk_shares.len(), n);
    }
}
