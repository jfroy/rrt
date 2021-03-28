use super::camera::*;
use super::materials::*;
use super::scene::*;
use super::sphere::*;
use super::types::*;
use crate::acceleration::BvhPartitionMethod;

pub fn chap11_scene<'a>(nx: usize, ny: usize) -> (Scene, Camera) {
    let mut scene = Scene::new();

    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(0., 0., -1.),
        radius: 0.5,
        material: Box::new(Lambertian {
            albedo: Vec4f::new(0.1, 0.2, 0.5, 1.),
        }),
    }));
    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(0., -100.5, -1.),
        radius: 100.,
        material: Box::new(Lambertian {
            albedo: Vec4f::new(0.8, 0.8, 0., 1.),
        }),
    }));
    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(1., 0., -1.),
        radius: 0.5,
        material: Box::new(Metal {
            albedo: Vec4f::new(0.8, 0.6, 0.2, 1.),
            fuzz: 0.3,
        }),
    }));
    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(-1., 0., -1.),
        radius: 0.5,
        material: Box::new(Dielectric { ref_idx: 1.5 }),
    }));
    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(-1., 0., -1.),
        radius: -0.45,
        material: Box::new(Dielectric { ref_idx: 1.5 }),
    }));

    scene.build_bvh(BvhPartitionMethod::Middle);

    let look_from = Vec4f::new_direction(3., 3., 2.);
    let look_at = Vec4f::new_direction(0., 0., -1.);
    let focus_dist = (look_from - look_at).magnitude();
    let aperture = 2.;
    let camera = Camera::from(CameraCreateInfo {
        look_from,
        look_at,
        up: Vec4f::new_direction(0., 1., 0.),
        fov: 20.,
        aspect: nx as f32 / ny as f32,
        aperture,
        focus_dist,
        time0: 0.,
        time1: 0.,
    });

    (scene, camera)
}
