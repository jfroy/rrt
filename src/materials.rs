use super::rng::*;
use super::types::*;

pub struct Scattered {
  pub r: Ray,
  pub attenuation: Vec4f,
}

pub trait Material {
  fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<Scattered>;
}

pub struct Hit<'scene> {
  pub t: f32,
  pub p: Vec4f,
  pub normal: Vec4f,
  pub material: &'scene (dyn Material + Sync),
}

// Lambertian

pub struct Lambertian {
  pub albedo: Vec4f,
}

impl Material for Lambertian {
  fn scatter(&self, _: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<Scattered> {
    let target = hit.p + hit.normal + Vec4f::in_unit_sphere(rng);
    let origin = hit.p;
    let direction = target - hit.p;
    let r = Ray { origin, direction };
    let attenuation = self.albedo;
    Some(Scattered { r, attenuation })
  }
}

// Metal

pub struct Metal {
  pub albedo: Vec4f,
  pub fuzz: f32,
}

impl Material for Metal {
  fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<Scattered> {
    let origin = hit.p;
    let direction =
      r.direction.normalized().reflected(hit.normal) + self.fuzz * Vec4f::in_unit_sphere(rng);
    let r = Ray { origin, direction };
    let attenuation = self.albedo;
    if direction.dot(hit.normal) > 0. {
      return Some(Scattered { r, attenuation });
    }
    None
  }
}

// Dielectric

fn schlick(cosine: f32, ref_idx: f32) -> f32 {
  let r0 = ((1. - ref_idx) / (1. + ref_idx)).powf(2.);
  r0 + (1. - r0) * (1. - cosine).powf(5.)
}

pub struct Dielectric {
  pub ref_idx: f32,
}

impl Material for Dielectric {
  fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<Scattered> {
    let unit_direction = r.direction.normalized();
    let reflected = unit_direction.reflected(hit.normal);
    let dir_dot_normal = r.direction.dot(hit.normal);
    let (outward_normal, ni_over_nt, cosine) = if dir_dot_normal > 0. {
      (
        -hit.normal,
        self.ref_idx,
        self.ref_idx * dir_dot_normal / r.direction.magnitude(),
      )
    } else {
      (
        hit.normal,
        1. / self.ref_idx,
        -dir_dot_normal / r.direction.magnitude(),
      )
    };
    let refracted = unit_direction.refracted(outward_normal, ni_over_nt);
    let reflect_prob = if refracted != Vec4f::zero() {
      schlick(cosine, self.ref_idx)
    } else {
      1.
    };
    let r = if rng.gen::<f32>() < reflect_prob {
      Ray {
        origin: hit.p,
        direction: reflected,
      }
    } else {
      Ray {
        origin: hit.p,
        direction: refracted,
      }
    };
    let attenuation = Vec4f::new(1., 1., 1., 0.);
    Some(Scattered { r, attenuation })
  }
}
