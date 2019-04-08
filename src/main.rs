extern crate image;
extern crate nalgebra as na;
extern crate nalgebra_glm as glm;
extern crate palette;
extern crate rand;

use palette::{LinSrgb, Pixel, Srgb};
use rand::distributions::Uniform;
use rand::prelude::*;
use std::cell::RefCell;
use std::option::Option;

struct Hit {
  t: f32,
  p: glm::Vec3,
  normal: glm::Vec3,
}

impl Default for Hit {
  fn default() -> Self {
    Hit {
      t: 0f32,
      p: glm::Vec3::zeros(),
      normal: glm::Vec3::zeros(),
    }
  }
}

trait Hittable {
  fn hit(&self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit>;
}

struct Sphere {
  center: glm::Vec3,
  radius: f32,
}

impl Hittable for Sphere {
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
      });
    }
    let root = root_term_1 + root_term_2;
    if root < t_max && root > t_min {
      let hit_p = r.point_at_parameter(root);
      return Some(Hit {
        t: root,
        p: hit_p,
        normal: (hit_p - self.center) / self.radius,
      });
    }
    None
  }
}

struct Scene {
  spheres: Vec<Sphere>,
}

impl Hittable for Scene {
  fn hit(&self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
    let mut closest_t = t_max;
    let mut closest_hit: Option<Hit> = None;
    for sphere in self.spheres.iter() {
      if let Some(hit) = sphere.hit(r, t_min, closest_t) {
        closest_t = hit.t;
        closest_hit = Some(hit);
      }
    }
    closest_hit
  }
}

struct Ray {
  origin: glm::Vec3,
  direction: glm::Vec3,
}

impl Ray {
  fn point_at_parameter(&self, t: f32) -> glm::Vec3 {
    self.origin + (t * self.direction)
  }
}

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

fn color(r: &Ray, scene: &Scene) -> glm::Vec3 {
  if let Some(hit) = scene.hit(r, 0.001f32, std::f32::MAX) {

    // Scale and bias the normal from [-1, 1] to [0, 1] and interpret as sRGB.
    // return (hit.normal + glm::Vec3::repeat(1f32)) * 0.5f32;

    // Diffuse material.
    let target = hit.p + hit.normal + random_in_unit_sphere();
    return 0.5f32
      * color(
        &Ray {
          origin: hit.p,
          direction: target - hit.p,
        },
        scene,
      );
  }

  let white: glm::Vec3 = glm::Vec3::repeat(1f32);
  let sky_blue: glm::Vec3 = glm::vec3(0.5f32, 0.7f32, 1.0f32);
  let unit_direction = r.direction.normalize();
  let t = 0.5f32 * (unit_direction.y + 1.0f32);
  return ((1.0f32 - t) * white) + (t * sky_blue);
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
  let dim = glm::vec2(200f32, 100f32);
  let ns = 100;
  let mut imgbuf = image::ImageBuffer::new(dim.x as u32, dim.y as u32);

  let mut scene: Scene = Scene {
    spheres: Vec::new(),
  };
  scene.spheres.push(Sphere {
    center: glm::vec3(0f32, 0f32, -1f32),
    radius: 0.5f32,
  });
  scene.spheres.push(Sphere {
    center: glm::vec3(0f32, -100.5f32, -1f32),
    radius: 100f32,
  });

  let cam = Camera {
    lower_left_corner: glm::vec3(-2f32, -1f32, -1f32),
    horizontal: glm::vec3(4f32, 0f32, 0f32),
    vertical: glm::vec3(0f32, 2f32, 0f32),
    origin: glm::vec3(0f32, 0f32, 0f32),
  };

  let mut rng = SmallRng::from_entropy();

  for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
    let mut c = glm::Vec3::zeros();
    for _ in 0..ns {
      let uv = glm::vec2(
        x as f32 + rng.gen::<f32>(),
        dim.y - (y + 1) as f32 + rng.gen::<f32>(),
      )
      .component_div(&dim);
      let r = cam.gen_ray(uv);
      c += color(&r, &scene);
    }
    // The book uses a simple gamma 2.0 function, not the sRGB OETF.
    c.apply(|e| (e / (ns as f32)).sqrt() * 255.99f32);
    *pixel = image::Rgb([c.x as i32 as u8, c.y as i32 as u8, c.z as i32 as u8]);

    // NOTE: This code outputs proper sRGB.
    // c = c.component_div(&glm::Vec3::repeat(ns as f32));
    // let linc = LinSrgb::new(c.x, c.y, c.z);
    // let srgbc: [u8; 3] = Srgb::from_linear(linc).into_format().into_raw();
    // *pixel = image::Rgb(srgbc);
  }

  imgbuf.save("o.ppm").unwrap();
}
