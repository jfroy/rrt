#![allow(clippy::just_underscores_and_digits)]

extern crate image;
extern crate rand;
extern crate rayon;

use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::borrow::Borrow;
use std::cell::UnsafeCell;
use std::option::Option;
use std::ptr::NonNull;
use std::sync::Mutex;
use vek::vec::repr_simd::*;
use vek::Lerp;

type Vec3f = vec3::Vec3<f32>;
type Vec4f = vec4::Vec4<f32>;
type Rgbf32 = rgb::Rgb<f32>;

// RNG

/// The [xoshiro](http://xoshiro.di.unimi.it/) generator is particularly well
/// suited for path tracing: it is a best-in-class PRNG (from a statistical and
/// performance POV), it supports an efficient jump-ahead operation which is
/// essential to prevent threads from having similar patterns, and its
/// implementation supports simd generation (with rust nightly) which is useful
/// to generate random vectors.
pub type RttRng = rand_xoshiro::Xoshiro128Plus;

// vek::Vec traits

trait RngVector {
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
  #[inline]
  fn new_rng_direction(rng: &mut RttRng) -> Vec4f {
    Vec4f::new(rng.gen(), rng.gen(), rng.gen(), 0.)
  }

  #[inline]
  fn in_unit_sphere(rng: &mut RttRng) -> Vec4f {
    loop {
      let v = 2. * Vec4f::new_rng_direction(rng) - Vec4f::new_direction(1., 1., 1.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }

  #[inline]
  fn in_unit_disc(rng: &mut RttRng) -> Vec4f {
    loop {
      let v = 2. * Vec4f::new(rng.gen(), rng.gen(), 0., 0.) - Vec4f::new(1., 1., 0., 0.);
      if v.dot(v) < 1. {
        return v;
      }
    }
  }
}

// Ray

struct Ray {
  origin: Vec4f,
  direction: Vec4f,
}

impl Ray {
  fn point_at_parameter(&self, t: f32) -> Vec4f {
    self.origin + (t * self.direction)
  }
}

// Material

struct Scattered {
  r: Ray,
  attenuation: Vec4f,
}

trait Material {
  fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<Scattered>;
}

// Hittable

struct Hit<'scene> {
  t: f32,
  p: Vec4f,
  normal: Vec4f,
  material: &'scene (dyn Material + Sync),
}

trait Hittable {
  fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>>;
}

// Sphere

struct Sphere {
  center: Vec4f,
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
  albedo: Vec4f,
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

struct Metal {
  albedo: Vec4f,
  fuzz: f32,
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

#[inline]
fn schlick(cosine: f32, ref_idx: f32) -> f32 {
  let r0 = ((1. - ref_idx) / (1. + ref_idx)).powf(2.);
  r0 + (1. - r0) * (1. - cosine).powf(5.)
}

struct Dielectric {
  ref_idx: f32,
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

// Camera

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
  fn new(
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

  fn gen_ray(&self, s: f32, t: f32, rng: &mut RttRng) -> Ray {
    let rd = self.lens_radius * Vec4f::in_unit_disc(rng);
    let offset = (self.u * rd.x) + (self.v * rd.y);
    let origin = self.origin + offset;
    Ray {
      origin,
      direction: self.lower_left_corner + (s * self.horizontal) + (t * self.vertical) - origin,
    }
  }
}

// Scene

pub struct Scene {
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

// Traces a ray. This is `color` in the book.
fn trace(r: &Ray, scene: &Scene, depth: i32, rng: &mut RttRng) -> Rgbf32 {
  if let Some(hit) = scene.hit(r, 0.001, std::f32::MAX) {
    if depth >= 50 {
      return Rgbf32::black();
    }
    if let Some(sc) = hit.material.scatter(r, &hit, rng) {
      return Rgbf32::from(Vec3f::from(sc.attenuation)) * trace(&sc.r, scene, depth + 1, rng);
    }
    return Rgbf32::black();
  }
  const WHITE: Rgbf32 = Rgbf32::new(1., 1., 1.);
  const SKY_BLUE: Rgbf32 = Rgbf32::new(0.5, 0.7, 1.);
  let unit_direction = r.direction.normalized();
  let t = 0.5 * (unit_direction.y + 1.);
  Lerp::lerp(WHITE, SKY_BLUE, t)
}

#[allow(dead_code)]
pub fn chap11_scene(nx: usize, ny: usize) -> (Scene, Camera) {
  let mut scene = Scene { spheres: vec![] };

  scene.spheres.push(Sphere {
    center: Vec4f::new_point(0., 0., -1.),
    radius: 0.5,
    material: Box::new(Lambertian {
      albedo: Vec4f::new(0.1, 0.2, 0.5, 1.),
    }),
  });
  scene.spheres.push(Sphere {
    center: Vec4f::new_point(0., -100.5, -1.),
    radius: 100.,
    material: Box::new(Lambertian {
      albedo: Vec4f::new(0.8, 0.8, 0., 1.),
    }),
  });
  scene.spheres.push(Sphere {
    center: Vec4f::new_point(1., 0., -1.),
    radius: 0.5,
    material: Box::new(Metal {
      albedo: Vec4f::new(0.8, 0.6, 0.2, 1.),
      fuzz: 0.3,
    }),
  });
  scene.spheres.push(Sphere {
    center: Vec4f::new_point(-1., 0., -1.),
    radius: 0.5,
    material: Box::new(Dielectric { ref_idx: 1.5 }),
  });
  scene.spheres.push(Sphere {
    center: Vec4f::new_point(-1., 0., -1.),
    radius: -0.45,
    material: Box::new(Dielectric { ref_idx: 1.5 }),
  });

  let look_from = Vec4f::new_point(3., 3., 2.);
  let look_at = Vec4f::new_point(0., 0., -1.);
  let dist_to_focus = (look_from - look_at).magnitude();
  let aperture = 2.;
  let camera = Camera::new(
    look_from,
    look_at,
    Vec4f::new_direction(0., 1., 0.),
    20.,
    nx as f32 / ny as f32,
    aperture,
    dist_to_focus,
  );

  (scene, camera)
}

pub fn chap12_scene(nx: usize, ny: usize, rng: &mut RttRng) -> (Scene, Camera) {
  let mut scene = Scene { spheres: vec![] };

  scene.spheres.push(Sphere {
    center: Vec4f::new_point(0., -1000., 0.),
    radius: 1000.,
    material: Box::new(Lambertian {
      albedo: Vec4f::new(0.5, 0.5, 0.5, 1.),
    }),
  });

  for a in -11..11 {
    for b in -11..11 {
      let center = Vec4f::new_point(
        a as f32 + 0.9 * rng.gen::<f32>(),
        0.2,
        b as f32 + 0.9 * rng.gen::<f32>(),
      );
      if (center - Vec4f::new_point(4., 0.2, 0.)).magnitude() > 0.9 {
        let choose_mat = rng.gen::<f32>();
        let sphere = if choose_mat < 0.8 {
          // Diffuse
          Sphere {
            center,
            radius: 0.2,
            material: Box::new(Lambertian {
              albedo: Vec4f::new(
                rng.gen::<f32>() * rng.gen::<f32>(),
                rng.gen::<f32>() * rng.gen::<f32>(),
                rng.gen::<f32>() * rng.gen::<f32>(),
                1.,
              ),
            }),
          }
        } else if choose_mat < 0.95 {
          // Metal
          Sphere {
            center,
            radius: 0.2,
            material: Box::new(Metal {
              albedo: Vec4f::new(
                0.5 * (1. + rng.gen::<f32>()),
                0.5 * (1. + rng.gen::<f32>()),
                0.5 * (1. + rng.gen::<f32>()),
                1.,
              ),
              fuzz: 0.5 * rng.gen::<f32>(),
            }),
          }
        } else {
          // Glass
          Sphere {
            center,
            radius: 0.2,
            material: Box::new(Dielectric { ref_idx: 1.5 }),
          }
        };
        scene.spheres.push(sphere);
      }
    }
  }

  scene.spheres.push(Sphere {
    center: Vec4f::new_point(0., 1., 0.),
    radius: 1.,
    material: Box::new(Dielectric { ref_idx: 1.5 }),
  });
  scene.spheres.push(Sphere {
    center: Vec4f::new_point(-4., 1., 0.),
    radius: 1.,
    material: Box::new(Lambertian {
      albedo: Vec4f::new_point(0.4, 0.2, 0.1),
    }),
  });
  scene.spheres.push(Sphere {
    center: Vec4f::new_point(4., 1., 0.),
    radius: 1.,
    material: Box::new(Metal {
      albedo: Vec4f::new_point(0.7, 0.6, 0.5),
      fuzz: 0.,
    }),
  });

  let look_from = Vec4f::new_point(13., 2., 3.);
  let look_at = Vec4f::new_point(0., 0., 0.);
  let up = Vec4f::new_direction(0., 1., 0.);
  let fov = 20.;
  let aspect = nx as f32 / ny as f32;
  let aperture = 0.1;
  let focus_dist = 10.;
  let camera = Camera::new(look_from, look_at, up, fov, aspect, aperture, focus_dist);

  (scene, camera)
}

thread_local!(static THREAD_RNG_KEY: UnsafeCell<RttRng> = UnsafeCell::new(RttRng::seed_from_u64(0)));

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

pub fn tracescene(
  nx: usize,
  ny: usize,
  ns: usize,
  scene: &Scene,
  camera: &Camera,
  pool: &rayon::ThreadPool,
) -> Vec<u8> {
  const BYTES_PER_PIXEL: usize = 3;
  let mut pixels = vec![0u8; ny * nx * BYTES_PER_PIXEL];
  pool.install(|| {
    pixels
      .par_chunks_mut(BYTES_PER_PIXEL)
      .enumerate()
      .for_each(|(idx, chunk)| {
        let raw = THREAD_RNG_KEY.with(|uc| uc.get());
        let mut nn = NonNull::new(raw).unwrap();
        let mut rng = unsafe { nn.as_mut() };
        let x = (idx % nx) as f32;
        let y = (ny - 1 - idx / nx) as f32;
        let mut c = Rgbf32::black();
        for _ in 0..ns {
          let ray = camera.gen_ray(
            (x + rng.gen::<f32>()) / nx as f32,
            (y + rng.gen::<f32>()) / ny as f32,
            &mut rng,
          );
          c += trace(&ray, scene, 0, &mut rng);
        }
        // The book uses a simple gamma 2.0 function, not the sRGB OETF.
        c.apply(|e| (e / (ns as f32)).sqrt() * 255.99);
        chunk[0] = c.r as u8;
        chunk[1] = c.g as u8;
        chunk[2] = c.b as u8;
      });
  });
  pixels
}
