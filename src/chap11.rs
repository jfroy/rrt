use super::camera::*;
use super::materials::*;
use super::scene::*;
use super::types::*;

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
