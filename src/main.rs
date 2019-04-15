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
    if discriminant < 0f32 {
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

fn trace(r: &Ray, scene: &Scene, depth: i32) -> glm::Vec3 {
  if let Some(hit) = scene.hit(r, 0.001f32, std::f32::MAX) {
    if depth >= 50 {
      return glm::Vec3::zeros();
    }
    if let Some(sc) = hit.material.scatter(r, &hit) {
      return glm::matrix_comp_mult(&sc.attenuation, &trace(&sc.r, scene, depth + 1));
    }
    return glm::Vec3::zeros();
  }

  let white: glm::Vec3 = glm::Vec3::repeat(1f32);
  let sky_blue: glm::Vec3 = glm::vec3(0.5f32, 0.7f32, 1.0f32);
  let unit_direction = r.direction.normalize();
  let t = 0.5f32 * (unit_direction.y + 1.0f32);
  return glm::lerp(&white, &sky_blue, t);
}

thread_local!(static RNG: RefCell<SmallRng> = RefCell::new(SmallRng::from_entropy()));

#[inline(always)]
fn random_in_unit_sphere() -> glm::Vec3 {
  let d = Uniform::from(0f32..1f32);
  loop {
    let v = RNG.with(|rng_rc| glm::Vec3::from_distribution(&d, &mut *rng_rc.borrow_mut()));
    let p = 2.0f32 * v - glm::Vec3::repeat(1f32);
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

  let lamb1 = Lambertian {
    albedo: glm::vec3(0.8f32, 0.3f32, 0.3f32),
  };
  let lamb2 = Lambertian {
    albedo: glm::vec3(0.8f32, 0.8f32, 0f32),
  };

  scene.spheres.push(Sphere {
    center: glm::vec3(0f32, 0f32, -1f32),
    radius: 0.5f32,
    material: &lamb1,
  });
  scene.spheres.push(Sphere {
    center: glm::vec3(0f32, -100.5f32, -1f32),
    radius: 100f32,
    material: &lamb2,
  });

  let cam = Camera {
    lower_left_corner: glm::vec3(-2f32, -1f32, -1f32),
    horizontal: glm::vec3(4f32, 0f32, 0f32),
    vertical: glm::vec3(0f32, 2f32, 0f32),
    origin: glm::vec3(0f32, 0f32, 0f32),
  };

  let rd = Uniform::from(0f32..1f32);

  let mut pixcoords: Vec<(f32, f32)> = vec![];
  for y in 0..dim.y {
    for x in 0..dim.x {
      pixcoords.push((y as f32, x as f32));
    }
  }
  let pixels: Vec<u8> = pixcoords
    .par_iter()
    .flat_map(|(y, x)| {
      let mut c = glm::Vec3::zeros();
      for _ in 0..ns {
        let rv = RNG.with(|rng_rc| glm::Vec2::from_distribution(&rd, &mut *rng_rc.borrow_mut()));
        let uv = (glm::vec2(*x, dim.y as f32 - *y + 1f32) + rv)
          .component_div(&glm::vec2(dim.x as f32, dim.y as f32));
        let r = cam.gen_ray(uv);
        c += trace(&r, &scene, 0);
      }
      // The book uses a simple gamma 2.0 function, not the sRGB OETF.
      c.apply(|e| (e / (ns as f32)).sqrt() * 255.99f32);
      vec![c.x as i32 as u8, c.y as i32 as u8, c.z as i32 as u8]
    })
    .collect();

  image::save_buffer(
    "o.ppm",
    &pixels[..],
    dim.x as u32,
    dim.y as u32,
    image::RGB(8),
  )
  .unwrap()
}
