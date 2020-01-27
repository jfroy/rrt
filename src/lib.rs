#![warn(clippy::all)]

extern crate image;
extern crate palette;
extern crate rand;
extern crate rayon;

use rand::prelude::*;
use rayon::prelude::*;
use std::borrow::Borrow;
use std::cell::RefCell;
use std::option::Option;
use vek::Lerp;

type Vec3f = vek::vec::Vec3<f32>;

// RNG

thread_local!(static RNG: RefCell<SmallRng> = RefCell::new(SmallRng::from_entropy()));

// vek::Vec traits

trait RngVector {
  /// Makes a new random vector using `rng`.
  ///
  /// Each component will be in the half-open range `[0,1)` (see
  /// https://rust-num.github.io/num/rand/trait.Rng.html#method.next_f32).
  ///
  /// This function does *not* create unit vectors.
  fn with_rng(rng: &mut impl Rng) -> Self;

  /// Generates a random vector inside the unit sphere.
  fn in_unit_sphere(rng: &mut impl Rng) -> Self;

  /// Generates a random vector inside the unit disc in the XY plane. The Z
  /// component shall be 0.
  fn in_unit_disc(rng: &mut impl Rng) -> Self;
}

impl RngVector for Vec3f {
  #[inline]
  fn with_rng(rng: &mut impl Rng) -> Vec3f {
    Vec3f::new(rng.gen(), rng.gen(), rng.gen())
  }

  #[inline]
  fn in_unit_sphere(rng: &mut impl Rng) -> Vec3f {
    loop {
      let v = 2. * Vec3f::with_rng(rng) - Vec3f::broadcast(1.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }

  #[inline]
  fn in_unit_disc(rng: &mut impl Rng) -> Vec3f {
    loop {
      let v = 2. * Vec3f::new(rng.gen(), rng.gen(), 0.) - Vec3f::new(1., 1., 0.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }
}

// Ray

struct Ray {
  origin: Vec3f,
  direction: Vec3f,
}

impl Ray {
  fn point_at_parameter(&self, t: f32) -> Vec3f {
    self.origin + (t * self.direction)
  }
}

// Material

struct Scattered {
  r: Ray,
  attenuation: Vec3f,
}

trait Material {
  fn scatter(&self, r: &Ray, hit: &Hit) -> Option<Scattered>;
}

// Hittable

struct Hit<'scene> {
  t: f32,
  p: Vec3f,
  normal: Vec3f,
  material: &'scene (dyn Material + Sync),
}

trait Hittable {
  fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>>;
}

// Sphere

struct Sphere {
  center: Vec3f,
  radius: f32,
  material: Box<dyn Material + Sync>,
}

impl Hittable for Sphere {
  fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>> {
    let oc = r.origin - self.center;
    let a = r.direction.dot(r.direction);
    let b = oc.dot(r.direction);
    let c = oc.dot(oc) - (self.radius * self.radius);
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
        material: self.material.borrow(),
      });
    }
    let root = root_term_1 + root_term_2;
    if root < t_max && root > t_min {
      let hit_p = r.point_at_parameter(root);
      return Some(Hit {
        t: root,
        p: hit_p,
        normal: (hit_p - self.center) / self.radius,
        material: self.material.borrow(),
      });
    }
    None
  }
}

// Lambertian

struct Lambertian {
  albedo: Vec3f,
}

impl Material for Lambertian {
  fn scatter(&self, _: &Ray, hit: &Hit) -> Option<Scattered> {
    let target =
      hit.p + hit.normal + RNG.with(|rng_rc| Vec3f::in_unit_sphere(&mut *rng_rc.borrow_mut()));
    let origin = hit.p;
    let direction = target - hit.p;
    let r = Ray { origin, direction };
    let attenuation = self.albedo;
    Some(Scattered { r, attenuation })
  }
}

// Metal

struct Metal {
  albedo: Vec3f,
  fuzz: f32,
}

impl Material for Metal {
  fn scatter(&self, r: &Ray, hit: &Hit) -> Option<Scattered> {
    let origin = hit.p;
    let direction = r.direction.normalized().reflected(hit.normal)
      + self.fuzz * RNG.with(|rng_rc| Vec3f::in_unit_sphere(&mut *rng_rc.borrow_mut()));
    let r = Ray { origin, direction };
    let attenuation = self.albedo;
    if direction.dot(hit.normal) > 0. {
      return Some(Scattered { r, attenuation });
    }
    None
  }
}

// Dielectric

#[inline]
fn schlick(cosine: f32, ref_idx: f32) -> f32 {
  let r0 = ((1. - ref_idx) / (1. + ref_idx)).powf(2.);
  r0 + (1. - r0) * (1. - cosine).powf(5.)
}

struct Dielectric {
  ref_idx: f32,
}

impl Material for Dielectric {
  fn scatter(&self, r: &Ray, hit: &Hit) -> Option<Scattered> {
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
    let reflect_prob = if refracted != Vec3f::zero() {
      schlick(cosine, self.ref_idx)
    } else {
      1.
    };
    let r = if RNG.with(|rng_rc| rng_rc.borrow_mut().gen::<f32>()) < reflect_prob {
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
    let attenuation = Vec3f::broadcast(1.);
    Some(Scattered { r, attenuation })
  }
}

// Camera

struct Camera {
  lower_left_corner: Vec3f,
  horizontal: Vec3f,
  vertical: Vec3f,
  origin: Vec3f,
  u: Vec3f,
  v: Vec3f,
  w: Vec3f,
  lens_radius: f32,
}

impl Camera {
  fn new(
    origin: Vec3f,
    look_at: Vec3f,
    up: Vec3f,
    fov: f32,
    aspect: f32,
    aperture: f32,
    focus_dist: f32,
  ) -> Camera {
    let theta = fov * std::f32::consts::PI / 180.;
    let half_height = f32::tan(theta / 2.);
    let half_width = aspect * half_height;
    let w = (origin - look_at).normalized();
    let u = up.cross(w).normalized();
    let v = w.cross(u);
    Camera {
      lower_left_corner: origin
        - (half_width * focus_dist * u)
        - (half_height * focus_dist * v)
        - focus_dist * w,
      horizontal: 2. * half_width * focus_dist * u,
      vertical: 2. * half_height * focus_dist * v,
      origin,
      u,
      v,
      w,
      lens_radius: aperture / 2.,
    }
  }

  fn gen_ray(&self, s: f32, t: f32) -> Ray {
    let rd = self.lens_radius * RNG.with(|rng_rc| Vec3f::in_unit_disc(&mut *rng_rc.borrow_mut()));
    let offset = (self.u * rd.x) + (self.v * rd.y);
    let origin = self.origin + offset;
    Ray {
      origin,
      direction: self.lower_left_corner + (s * self.horizontal) + (t * self.vertical) - origin,
    }
  }
}

// Scene

struct Scene {
  spheres: Vec<Sphere>,
}

impl Hittable for Scene {
  fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>> {
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

// main

// Traces a ray. This is `color` in the book.
fn trace(r: &Ray, scene: &Scene, depth: i32) -> Vec3f {
  if let Some(hit) = scene.hit(r, 0.001, std::f32::MAX) {
    if depth >= 50 {
      return Vec3f::zero();
    }
    if let Some(sc) = hit.material.scatter(r, &hit) {
      return sc.attenuation * trace(&sc.r, scene, depth + 1);
    }
    return Vec3f::zero();
  }
  const WHITE: Vec3f = Vec3f::new(1., 1., 1.);
  const SKY_BLUE: Vec3f = Vec3f::new(0.5, 0.7, 1.0);
  let unit_direction = r.direction.normalized();
  let t = 0.5 * (unit_direction.y + 1.);
  Lerp::lerp(WHITE, SKY_BLUE, t)
}

fn chap11_scene() -> Scene {
  let mut scene = Scene { spheres: vec![] };

  scene.spheres.push(Sphere {
    center: Vec3f::new(0., 0., -1.),
    radius: 0.5,
    material: Box::new(Lambertian {
      albedo: Vec3f::new(0.1, 0.2, 0.5),
    }),
  });
  scene.spheres.push(Sphere {
    center: Vec3f::new(0., -100.5, -1.),
    radius: 100.,
    material: Box::new(Lambertian {
      albedo: Vec3f::new(0.8, 0.8, 0.),
    }),
  });
  scene.spheres.push(Sphere {
    center: Vec3f::new(1., 0., -1.),
    radius: 0.5,
    material: Box::new(Metal {
      albedo: Vec3f::new(0.8, 0.6, 0.2),
      fuzz: 0.3,
    }),
  });
  scene.spheres.push(Sphere {
    center: Vec3f::new(-1., 0., -1.),
    radius: 0.5,
    material: Box::new(Dielectric { ref_idx: 1.5 }),
  });
  scene.spheres.push(Sphere {
    center: Vec3f::new(-1., 0., -1.),
    radius: -0.45,
    material: Box::new(Dielectric { ref_idx: 1.5 }),
  });

  scene
}

pub fn tracescene(nx: usize, ny: usize, ns: usize) -> Vec<u8> {
  let scene = chap11_scene();

  let look_from = Vec3f::new(3., 3., 2.);
  let look_at = Vec3f::new(0., 0., -1.);
  let dist_to_focus = (look_from - look_at).magnitude();
  let aperture = 2.;
  let cam = Camera::new(
    look_from,
    look_at,
    Vec3f::new(0., 1., 0.),
    20.,
    nx as f32 / ny as f32,
    aperture,
    dist_to_focus,
  );

  const BYTES_PER_PIXEL: usize = 3;
  let mut pixels = vec![0u8; ny * nx * BYTES_PER_PIXEL];
  pixels
    .par_chunks_mut(BYTES_PER_PIXEL)
    .enumerate()
    .for_each(|(idx, chunk)| {
      let x = (idx % nx) as f32;
      let y = (ny - 1 - idx / nx) as f32;
      let mut c = Vec3f::zero();
      for _ in 0..ns {
        let (rx, ry) = RNG.with(|rng_rc| {
          let mut rng = rng_rc.borrow_mut();
          (rng.gen::<f32>(), rng.gen::<f32>())
        });
        let r = cam.gen_ray((x + rx) / nx as f32, (y + ry) / ny as f32);
        c += trace(&r, &scene, 0);
      }
      // The book uses a simple gamma 2.0 function, not the sRGB OETF.
      c.apply(|e| (e / (ns as f32)).sqrt() * 255.99);
      chunk[0] = c.x as u8;
      chunk[1] = c.y as u8;
      chunk[2] = c.z as u8;
    });
  pixels
}
