use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_std::{One, UniformRand};
use batch_threshold::{dealer::Dealer, decryption::SecretKey, encryption::encrypt};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::thread_rng;

type E = Bls12_381;
type Fr = <E as Pairing>::ScalarField;
type G1 = <E as Pairing>::G1;

fn bench_encrypt(c: &mut Criterion) {
    let mut rng = thread_rng();

    let n = 1 << 4;
    let mut group = c.benchmark_group("encrypt");
    for size in 2..6 {
        // timing doesn't change since encryption is independent of batch_size
        // done as a sanity check

        let batch_size = 1 << size;

        let mut dealer = Dealer::<E>::new(batch_size, n, n / 2 - 1);
        let (crs, sk_shares) = dealer.setup(&mut rng);
        let pk = dealer.get_pk();

        let mut secret_key: Vec<SecretKey<E>> = Vec::new();
        for i in 0..n {
            secret_key.push(SecretKey::new(sk_shares[i]));
        }

        let msg = [1u8; 32];

        let hid = G1::rand(&mut rng);

        group.bench_with_input(
            BenchmarkId::from_parameter(batch_size),
            &(msg, hid, crs.htau, pk),
            |b, &inp| {
                b.iter(|| encrypt::<E>(inp.0, Fr::one(), inp.1, inp.2, inp.3));
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_encrypt);
criterion_main!(benches);
