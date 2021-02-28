use super::rng::*;
use std::cell::UnsafeCell;
use std::ptr::NonNull;
use std::sync::Mutex;

thread_local!(pub static THREAD_RNG_KEY: UnsafeCell<RttRng> = UnsafeCell::new(RttRng::seed_from_u64(0)));

// Create a rayon thread pool with a start handler that installs a suitable rng
// for the thread. Each thread's RNG is built from the provided RNG doing
// `thread_index` jumps.
pub fn init_pool_with_rng(rng: RttRng) -> rayon::ThreadPool {
  let rng_mutex = Mutex::new(rng);
  rayon::ThreadPoolBuilder::new()
    .start_handler(move |idx| {
      let raw = THREAD_RNG_KEY.with(|uc| uc.get());
      let mut nn = NonNull::new(raw).unwrap();
      let mut rng = rng_mutex.lock().unwrap().clone();
      for _ in 0..idx {
        rng.jump();
      }
      unsafe {
        *nn.as_mut() = rng;
      }
    })
    .build()
    .unwrap()
}
