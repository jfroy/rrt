// Provides random number generation.

use super::types::*;
use rand::distributions::Uniform;
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
  /// Generates a random vector inside the unit sphere from a uniform
  /// distribution.
  fn gen_uniform_random_in_unit_sphere(rng: &mut RttRng) -> Self;

  /// Generates a random vector inside the unit disc in the XY plane from a
  /// uniform distribution. The Z component shall be 0.
  fn gen_uniform_random_in_unit_disc(rng: &mut RttRng) -> Self;

  /// Generates a random unit vector from a uniform distribution.
  fn gen_uniform_random_unit(rng: &mut RttRng) -> Self;
}

impl RngVector for Vec4f {
  fn gen_uniform_random_in_unit_sphere(rng: &mut RttRng) -> Vec4f {
    let d = Uniform::new_inclusive(-1., 1.);
    loop {
      let v = Vec4f::new(rng.sample(d), rng.sample(d), rng.sample(d), 0.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }

  fn gen_uniform_random_in_unit_disc(rng: &mut RttRng) -> Vec4f {
    let d = Uniform::new_inclusive(-1., 1.);
    loop {
      let v = Vec4f::new(rng.sample(d), rng.sample(d), 0., 0.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }

  fn gen_uniform_random_unit(rng: &mut RttRng) -> Vec4f {
    let v = Vec4f::gen_uniform_random_in_unit_sphere(rng);
    v.normalized()
  }
}
