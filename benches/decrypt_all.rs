use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::UniformRand;
use batch_threshold::{
    dealer::Dealer,
    decryption::{decrypt_all, SecretKey},
    encryption::{encrypt, Ciphertext},
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::thread_rng;

type E = Bls12_381;
type Fr = <E as Pairing>::ScalarField;
type G1 = <E as Pairing>::G1;

fn bench_decrypt_all(c: &mut Criterion) {
    let mut rng = thread_rng();

    let n = 1 << 4;
    let mut group = c.benchmark_group("decrypt_all");
    group.sample_size(20);

    for size in 2..=10 {
        let batch_size = 1 << size;

        let mut dealer = Dealer::<E>::new(batch_size, n, n / 2 - 1);
        let (crs, sk_shares) = dealer.setup(&mut rng);
        let pk = dealer.get_pk();

        let mut secret_key: Vec<SecretKey<E>> = Vec::new();
        for i in 0..n {
            secret_key.push(SecretKey::new(sk_shares[i]));
        }
        let public_keys = secret_key.iter().map(|sk| sk.get_pk()).collect::<Vec<_>>();

        let msg = [1u8; 32];

        let hid = G1::rand(&mut rng);

        let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();

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

        // bench full decryption
        group.bench_with_input(
            BenchmarkId::from_parameter(batch_size),
            &(public_keys, partial_decryptions, ct, hid, crs),
            |b, inp| {
                b.iter(|| decrypt_all(&inp.0, &inp.1, &inp.2, inp.3, &inp.4));
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_decrypt_all);
criterion_main!(benches);
