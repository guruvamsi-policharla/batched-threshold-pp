use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::UniformRand;
use batch_threshold::{
    dealer::Dealer,
    decryption::{decrypt_all, SecretKey},
    encryption::{encrypt, Ciphertext},
};

type E = Bls12_381;
type Fr = <E as Pairing>::ScalarField;
type G1 = <E as Pairing>::G1;

fn main() {
    let mut rng = ark_std::test_rng();
    let batch_size = 1 << 5;
    let n = 1 << 3;

    let mut dealer = Dealer::<E>::new(batch_size, n, n / 2 - 1);
    let (crs, sk_shares) = dealer.setup(&mut rng);

    let mut secret_key: Vec<SecretKey<E>> = Vec::new();
    for i in 0..n {
        secret_key.push(SecretKey::new(sk_shares[i]));
    }
    let public_keys = secret_key.iter().map(|sk| sk.get_pk()).collect::<Vec<_>>();

    let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();

    let msg = [1u8; 32];
    let hid = G1::rand(&mut rng);
    let pk = dealer.get_pk();

    // generate ciphertexts for all points in tx_domain
    let mut ct: Vec<Ciphertext<E>> = Vec::new();
    for x in tx_domain.elements() {
        ct.push(encrypt::<E>(msg, x, hid, crs.htau, pk));
    }

    // generate partial decryptions
    let mut partial_decryptions: Vec<G1> = Vec::new();
    for i in 0..n {
        let partial_decryption = secret_key[i].partial_decrypt(&ct, hid, pk, &crs);
        partial_decryptions.push(partial_decryption);
    }

    let messages = decrypt_all(&public_keys, &partial_decryptions, &ct, hid, &crs);
    for i in 0..batch_size {
        assert_eq!(msg, messages[i]);
    }
}
