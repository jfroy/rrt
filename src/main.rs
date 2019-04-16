#![warn(clippy::all)]

extern crate image;
extern crate nalgebra as na;
extern crate nalgebra_glm as glm;
extern crate palette;
extern crate rand;
extern crate rayon;

use palette::{LinSrgb, Pixel, Srgb};
use rand::distributions::Uniform;
use rand::prelude::*;
use rayon::prelude::*;
use std::cell::RefCell;
use std::option::Option;

// RNG

thread_local!(static RNG: RefCell<SmallRng> = RefCell::new(SmallRng::from_entropy()));

// Ray

struct Ray {
  origin: glm::Vec3,
  direction: glm::Vec3,
}

impl Ray {
  fn point_at_parameter(&self, t: f32) -> glm::Vec3 {
    self.origin + (t * self.direction)
  }
}

// Material

struct Scattered {
  r: Ray,
  attenuation: glm::Vec3,
}

trait Material {
  fn scatter(&self, r: &Ray, hit: &Hit) -> Option<Scattered>;
}

// Hittable

struct Hit<'obj> {
  t: f32,
  p: glm::Vec3,
  normal: glm::Vec3,
  material: &'obj (Material + Sync),
}

trait Hittable<'obj> {
  fn hit(&self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit>;
}

// Sphere

struct Sphere<'obj> {
  center: glm::Vec3,
  radius: f32,
  material: &'obj (Material + Sync),
}

impl<'obj> Hittable<'obj> for Sphere<'obj> {
  fn hit(&self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
    let oc = r.origin - self.center;
    let a = glm::dot(&r.direction, &r.direction);
    let b = glm::dot(&oc, &r.direction);
    let c = glm::dot(&oc, &oc) - (self.radius * self.radius);
    let discriminant = (b * b) - (a * c);
    if discriminant < 0. {
      return None;
    }
    let root_term_1 = -b / a;
    let root_term_2 = (discriminant).sqrt() / a;
    let root = root_term_1 - root_term_2;
    if root < t_max && root > t_min {
      let hit_p = r.point_at_parameter(root);
      return Some(Hit {
        t: root,
        p: hit_p,
        normal: (hit_p - self.center) / self.radius,
        material: self.material,
      });
    }
    let root = root_term_1 + root_term_2;
    if root < t_max && root > t_min {
      let hit_p = r.point_at_parameter(root);
      return Some(Hit {
        t: root,
        p: hit_p,
        normal: (hit_p - self.center) / self.radius,
        material: self.material,
      });
    }
    None
  }
}

// Scene

struct Scene<'obj> {
  spheres: Vec<Sphere<'obj>>,
}

impl<'obj> Hittable<'obj> for Scene<'obj> {
  fn hit(&self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
    let mut closest_t = t_max;
    let mut closest_hit: Option<Hit> = None;
    for sphere in &self.spheres {
      if let Some(hit) = sphere.hit(r, t_min, closest_t) {
        closest_t = hit.t;
        closest_hit = Some(hit);
      }
    }
    closest_hit
  }
}

// Lambertian

struct Lambertian {
  albedo: glm::Vec3,
}

impl Material for Lambertian {
  fn scatter(&self, _: &Ray, hit: &Hit) -> Option<Scattered> {
    let target = hit.p + hit.normal + random_in_unit_sphere();
    let origin = hit.p;
    let direction = target - hit.p;
    let r = Ray { origin, direction };
    let attenuation = self.albedo;
    Some(Scattered { r, attenuation })
  }
}

// Metal

struct Metal {
  albedo: glm::Vec3,
  fuzz: f32,
}

impl Material for Metal {
  fn scatter(&self, r: &Ray, hit: &Hit) -> Option<Scattered> {
    let origin = hit.p;
    let direction =
      glm::reflect_vec(&r.direction.normalize(), &hit.normal) + self.fuzz * random_in_unit_sphere();
    let r = Ray { origin, direction };
    let attenuation = self.albedo;
    if glm::dot(&direction, &hit.normal) > 0. {
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

struct Dielectric {
  ref_idx: f32,
}

impl Material for Dielectric {
  fn scatter(&self, r: &Ray, hit: &Hit) -> Option<Scattered> {
    let unit_direction = r.direction.normalize();
    let reflected = glm::reflect_vec(&unit_direction, &hit.normal);
    let dir_dot_normal = glm::dot(&r.direction, &hit.normal);
    let (outward_normal, ni_over_nt, cosine) = if dir_dot_normal > 0. {
      (
        -hit.normal,
        self.ref_idx,
        self.ref_idx * dir_dot_normal / glm::length(&r.direction),
      )
    } else {
      (
        hit.normal,
        1. / self.ref_idx,
        -dir_dot_normal / glm::length(&r.direction),
      )
    };
    let refracted = glm::refract_vec(&unit_direction, &outward_normal, ni_over_nt);
    let reflect_prob = if refracted != glm::Vec3::zeros() {
      schlick(cosine, self.ref_idx)
    } else {
      1.
    };
    let r = if RNG.with(|rng_rc| Uniform::from(0f32..1f32).sample(&mut *rng_rc.borrow_mut()))
      < reflect_prob
    {
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
    let attenuation = glm::vec3(1., 1., 1.);
    Some(Scattered { r, attenuation })
  }
}

// Camera

struct Camera {
  lower_left_corner: glm::Vec3,
  horizontal: glm::Vec3,
  vertical: glm::Vec3,
  origin: glm::Vec3,
}

impl Camera {
  fn gen_ray(&self, uv: glm::Vec2) -> Ray {
    Ray {
      origin: self.origin,
      direction: self.lower_left_corner + (uv.x * self.horizontal) + (uv.y * self.vertical)
        - self.origin,
    }
  }
}

// main

// Traces a ray. This is `color` in the book.
fn trace(r: &Ray, scene: &Scene, depth: i32) -> glm::Vec3 {
  if let Some(hit) = scene.hit(r, 0.001, std::f32::MAX) {
    if depth >= 50 {
      return glm::Vec3::zeros();
    }
    if let Some(sc) = hit.material.scatter(r, &hit) {
      return glm::matrix_comp_mult(&sc.attenuation, &trace(&sc.r, scene, depth + 1));
    }
    return glm::Vec3::zeros();
  }
  let white: glm::Vec3 = glm::Vec3::repeat(1.);
  let sky_blue: glm::Vec3 = glm::vec3(0.5, 0.7, 1.0);
  let unit_direction = r.direction.normalize();
  let t = 0.5 * (unit_direction.y + 1.);
  glm::lerp(&white, &sky_blue, t)
}

#[inline(always)]
fn random_in_unit_sphere() -> glm::Vec3 {
  let d = Uniform::from(0f32..1f32);
  loop {
    let v = RNG.with(|rng_rc| glm::Vec3::from_distribution(&d, &mut *rng_rc.borrow_mut()));
    let p = 2. * v - glm::Vec3::repeat(1.);
    if glm::dot(&p, &p) < 1.0 {
      return p;
    }
  }

  // NOTE: This produces a much smaller shadow for the small sphere and other
  // differences from the reference implementation that generates vectors inside
  // the unit sphere (as opposed to on the surface).
  // let sphere = rand::distributions::UnitSphereSurface::new();
  // let s = RNG.with(|rng_rc| sphere.sample(&mut *rng_rc.borrow_mut()));
  // glm::vec3(s[0] as f32, s[1] as f32, s[2] as f32)
}

fn main() {
  let dim = glm::vec2(1000, 500);
  let ns = 100;

  let mut scene: Scene = Scene {
    spheres: Vec::new(),
  };

  let mat_lamb1 = Lambertian {
    albedo: glm::vec3(0.1, 0.2, 0.5),
  };
  let mat_lamb2 = Lambertian {
    albedo: glm::vec3(0.8, 0.8, 0.),
  };
  let mat_metal1 = Metal {
    albedo: glm::vec3(0.8, 0.6, 0.2),
    fuzz: 0.3,
  };
  let mat_dia1 = Dielectric { ref_idx: 1.5 };

  scene.spheres.push(Sphere {
    center: glm::vec3(0., 0., -1.),
    radius: 0.5,
    material: &mat_lamb1,
  });
  scene.spheres.push(Sphere {
    center: glm::vec3(0., -100.5, -1.),
    radius: 100.,
    material: &mat_lamb2,
  });
  scene.spheres.push(Sphere {
    center: glm::vec3(1., 0., -1.),
    radius: 0.5,
    material: &mat_metal1,
  });
  scene.spheres.push(Sphere {
    center: glm::vec3(-1., 0., -1.),
    radius: 0.5,
    material: &mat_dia1,
  });
  scene.spheres.push(Sphere {
    center: glm::vec3(-1., 0., -1.),
    radius: -0.45,
    material: &mat_dia1,
  });

  let cam = Camera {
    lower_left_corner: glm::vec3(-2., -1., -1.),
    horizontal: glm::vec3(4., 0., 0.),
    vertical: glm::vec3(0., 2., 0.),
    origin: glm::vec3(0., 0., 0.),
  };

  let rd = Uniform::from(0f32..1f32);

  const BYTES_PER_PIXEL: usize = 3;
  let mut pixels = vec![0u8; dim.y as usize * dim.x as usize * BYTES_PER_PIXEL];
  pixels
    .par_chunks_mut(BYTES_PER_PIXEL)
    .enumerate()
    .for_each(|(idx, chunk)| {
      let puv = glm::vec2((idx % dim.x) as f32, (dim.y - 1 - idx / dim.x) as f32);
      let mut c = glm::Vec3::zeros();
      for _ in 0..ns {
        let noise = RNG.with(|rng_rc| glm::Vec2::from_distribution(&rd, &mut *rng_rc.borrow_mut()));
        let suv = (puv + noise).component_div(&glm::vec2(dim.x as f32, dim.y as f32));
        let r = cam.gen_ray(suv);
        c += trace(&r, &scene, 0);
      }
      // The book uses a simple gamma 2.0 function, not the sRGB OETF.
      c.apply(|e| (e / (ns as f32)).sqrt() * 255.99);
      chunk[0] = c.x as u8;
      chunk[1] = c.y as u8;
      chunk[2] = c.z as u8;
    });

  image::save_buffer(
    "o.ppm",
    &pixels[..],
    dim.x as u32,
    dim.y as u32,
    image::RGB(8),
  )
  .unwrap()
}
