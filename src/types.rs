// Provides basic types.

use vek::vec::repr_simd::*;

pub type Vec3f = vec3::Vec3<f32>;
pub type Vec4f = vec4::Vec4<f32>;
pub type Rgbf32 = rgb::Rgb<f32>;

pub struct Ray {
  pub origin: Vec4f,
  pub direction: Vec4f,
}

impl Ray {
  pub fn point_at_parameter(&self, t: f32) -> Vec4f {
    self.origin + (t * self.direction)
  }
}
