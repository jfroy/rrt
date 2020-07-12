use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use rand::prelude::*;
use rrt::tracescene;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("tracescene/10x10x4", move |b| {
    const NX: usize = 10;
    const NY: usize = 10;
    const NS: usize = 4;
    let mut rng = rand::rngs::SmallRng::seed_from_u64(0xC0DE_D09E);
    // tracescene returns a large Vec (the image), so use iter_batched to deal with the
    // memory drop. iter_with_large_drop is not suitable because the outputs are too large
    // to accumulate.
    b.iter_batched(|| (), |_| tracescene(NX, NY, NS, &mut rng), BatchSize::SmallInput);
  });
}

criterion_group! {
  name = benches;
  config = Criterion::default().measurement_time(Duration::new(10, 0));
  targets = criterion_benchmark
}

criterion_main!(benches);
