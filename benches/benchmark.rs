use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use rand::SeedableRng;
use rrt::chap12::*;
use rrt::rng::*;
use rrt::threadpool::*;
use rrt::tracescene;
use std::sync::atomic::AtomicUsize;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
  const NX: usize = 10;
  const NY: usize = 10;
  const NS: usize = 4;
  let mut rng = RttRng::seed_from_u64(0);
  let (scene, camera) = chap12_scene(NX, NY, &mut rng);
  let pool = init_pool_with_rng(rng);
  c.bench_function("tracescene/10x10x4", move |b| {
    // tracescene returns a large Vec (the image), so use iter_batched to deal with the
    // memory drop. iter_with_large_drop is not suitable because the outputs are too large
    // to accumulate.
    b.iter_batched(
      || (),
      |_| tracescene(NX, NY, NS, &scene, &camera, &pool, &AtomicUsize::new(0)),
      BatchSize::SmallInput,
    );
  });
}

criterion_group! {
  name = benches;
  config = Criterion::default().measurement_time(Duration::new(10, 0));
  targets = criterion_benchmark
}

criterion_main!(benches);
