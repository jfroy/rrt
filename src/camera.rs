use super::rng::*;
use super::types::*;

pub struct Camera {
  lower_left_corner: Vec4f,
  horizontal: Vec4f,
  vertical: Vec4f,
  origin: Vec4f,
  u: Vec4f,
  v: Vec4f,
  #[allow(dead_code)]
  w: Vec4f,
  lens_radius: f32,
}

impl Camera {
  pub fn new(
    origin: Vec4f,
    look_at: Vec4f,
    up: Vec4f,
    fov: f32,
    aspect: f32,
    aperture: f32,
    focus_dist: f32,
  ) -> Camera {
    let theta = fov * std::f32::consts::PI / 180.;
    let half_height = f32::tan(theta / 2.);
    let half_width = aspect * half_height;
    let w = (origin - look_at).normalized();
    let u = Vec3f::from(up).cross(Vec3f::from(w)).normalized();
    let v = Vec3f::from(w).cross(u);
    Camera {
      lower_left_corner: origin
        - (half_width * focus_dist * u)
        - (half_height * focus_dist * v)
        - focus_dist * w,
      horizontal: Vec4f::from_direction(2. * half_width * focus_dist * u),
      vertical: Vec4f::from_direction(2. * half_height * focus_dist * v),
      origin,
      u: Vec4f::from_direction(u),
      v: Vec4f::from_direction(v),
      w,
      lens_radius: aperture / 2.,
    }
  }

  pub fn gen_ray(&self, s: f32, t: f32, rng: &mut RttRng) -> Ray {
    let rd = self.lens_radius * Vec4f::gen_uniform_random_in_unit_disc(rng);
    let offset = (self.u * rd.x) + (self.v * rd.y);
    let origin = self.origin + offset;
    Ray {
      origin,
      direction: self.lower_left_corner + (s * self.horizontal) + (t * self.vertical) - origin,
    }
  }
}
