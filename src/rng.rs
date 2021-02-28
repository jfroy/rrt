// Provides random number generation.

use super::types::*;
pub use rand::Rng;
pub use rand::SeedableRng;

/// The [xoshiro](http://xoshiro.di.unimi.it/) generator is particularly well
/// suited for path tracing: it is a best-in-class PRNG (from a statistical and
/// performance POV), it supports an efficient jump-ahead operation which is
/// essential to prevent threads from having similar patterns, and its
/// implementation supports simd generation (with rust nightly) which is useful
/// to generate random vectors.
pub type RttRng = rand_xoshiro::Xoshiro128Plus;

// vek::Vec traits

pub trait RngVector {
  /// Makes a new random direction vector using `rng`.
  ///
  /// Each component will be in the half-open range `[0,1)` (see
  /// https://rust-random.github.io/rand/rand/distributions/struct.Standard.html#floating-point-implementation.
  ///
  /// This function does *not* create unit vectors.
  fn new_rng_direction(rng: &mut RttRng) -> Self;

  /// Generates a random vector inside the unit sphere.
  fn in_unit_sphere(rng: &mut RttRng) -> Self;

  /// Generates a random vector inside the unit disc in the XY plane. The Z
  /// component shall be 0.
  fn in_unit_disc(rng: &mut RttRng) -> Self;
}

impl RngVector for Vec4f {
  fn new_rng_direction(rng: &mut RttRng) -> Vec4f {
    Vec4f::new(rng.gen(), rng.gen(), rng.gen(), 0.)
  }

  fn in_unit_sphere(rng: &mut RttRng) -> Vec4f {
    loop {
      let v = 2. * Vec4f::new_rng_direction(rng) - Vec4f::new_direction(1., 1., 1.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }

  fn in_unit_disc(rng: &mut RttRng) -> Vec4f {
    loop {
      let v = 2. * Vec4f::new(rng.gen(), rng.gen(), 0., 0.) - Vec4f::new(1., 1., 0., 0.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }
}
