use super::camera::*;
use super::materials::*;
use super::rng::*;
use super::scene::*;
use super::sphere::*;
use super::types::*;
use rand::distributions::Uniform;

pub fn book2_chap2_scene(nx: usize, ny: usize, rng: &mut RttRng) -> (Scene, Camera) {
    let mut scene = Scene::new();
    // Clone the RNG for the moving sphere center1 calculations to keep the scene
    // otherwise identical to chap12_scene.
    let mut rng2 = rng.clone();

    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(0., -1000., 0.),
        radius: 1000.,
        material: Box::new(Lambertian {
            albedo: Vec4f::new(0.5, 0.5, 0.5, 1.),
        }),
    }));

    for a in -11..11 {
        for b in -11..11 {
            let center = Vec4f::new_direction(
                a as f32 + 0.9 * rng.gen::<f32>(),
                0.2,
                b as f32 + 0.9 * rng.gen::<f32>(),
            );
            if (center - Vec4f::new_direction(4., 0.2, 0.)).magnitude() > 0.9 {
                let choose_mat = rng.gen::<f32>();
                let sphere = if choose_mat < 0.8 {
                    // Diffuse
                    let center0 = center;
                    let center1 = center0
                        + Vec4f::new_direction(
                            0.,
                            rng2.sample(Uniform::new_inclusive(0., 0.5)),
                            0.,
                        );
                    Sphere {
                        center0,
                        center1,
                        time0: 0.,
                        time1: 1.,
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
                    Sphere::from(StationarySphere {
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
                    })
                } else {
                    // Glass
                    Sphere::from(StationarySphere {
                        center,
                        radius: 0.2,
                        material: Box::new(Dielectric { ref_idx: 1.5 }),
                    })
                };
                scene.spheres.push(sphere);
            }
        }
    }

    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(0., 1., 0.),
        radius: 1.,
        material: Box::new(Dielectric { ref_idx: 1.5 }),
    }));
    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(-4., 1., 0.),
        radius: 1.,
        material: Box::new(Lambertian {
            albedo: Vec4f::new_direction(0.4, 0.2, 0.1),
        }),
    }));
    scene.spheres.push(Sphere::from(StationarySphere {
        center: Vec4f::new_direction(4., 1., 0.),
        radius: 1.,
        material: Box::new(Metal {
            albedo: Vec4f::new_direction(0.7, 0.6, 0.5),
            fuzz: 0.,
        }),
    }));

    let look_from = Vec4f::new_direction(13., 2., 3.);
    let look_at = Vec4f::new_direction(0., 0., 0.);
    let up = Vec4f::new_direction(0., 1., 0.);
    let fov = 20.;
    let aspect = nx as f32 / ny as f32;
    let aperture = 0.1;
    let focus_dist = 10.;
    let camera = Camera::from(CameraCreateInfo {
        look_from,
        look_at,
        up,
        fov,
        aspect,
        aperture,
        focus_dist,
        time0: 0.,
        time1: 1.,
    });

    (scene, camera)
}
