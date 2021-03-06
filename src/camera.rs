use super::rng::*;
use super::types::*;
use rand::distributions::Uniform;

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
  time_dist: Uniform<f32>,
}

pub struct CameraCreateInfo {
  pub look_from: Vec4f,
  pub look_at: Vec4f,
  pub up: Vec4f,
  pub fov: f32,
  pub aspect: f32,
  pub aperture: f32,
  pub focus_dist: f32,
  pub time0: f32,
  pub time1: f32,
}

impl Camera {
  pub fn from(ci: CameraCreateInfo) -> Camera {
    let theta = ci.fov * std::f32::consts::PI / 180.;
    let half_height = f32::tan(theta / 2.);
    let half_width = ci.aspect * half_height;
    let w = (ci.look_from - ci.look_at).normalized();
    let u = Vec3f::from(ci.up).cross(Vec3f::from(w)).normalized();
    let v = Vec3f::from(w).cross(u);
    Camera {
      lower_left_corner: ci.look_from
        - (half_width * ci.focus_dist * u)
        - (half_height * ci.focus_dist * v)
        - ci.focus_dist * w,
      horizontal: Vec4f::from_direction(2. * half_width * ci.focus_dist * u),
      vertical: Vec4f::from_direction(2. * half_height * ci.focus_dist * v),
      origin: ci.look_from,
      u: Vec4f::from_direction(u),
      v: Vec4f::from_direction(v),
      w,
      lens_radius: ci.aperture / 2.,
      time_dist: Uniform::new_inclusive(ci.time0, ci.time1),
    }
  }

  pub fn gen_ray(&self, s: f32, t: f32, rng: &mut RttRng) -> Ray {
    let rd = self.lens_radius * Vec4f::gen_uniform_random_in_unit_disc(rng);
    let offset = (self.u * rd.x) + (self.v * rd.y);
    let origin = self.origin + offset;
    Ray {
      origin,
      direction: self.lower_left_corner + (s * self.horizontal) + (t * self.vertical) - origin,
      time: rng.sample(self.time_dist),
    }
  }
}
