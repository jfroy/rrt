#![allow(clippy::just_underscores_and_digits)]
#![feature(const_fn_floating_point_arithmetic)]
#![feature(platform_intrinsics)]

pub mod book2chap2;
pub mod chap11;
pub mod chap12;
pub mod rng;
pub mod threadpool;

mod acceleration;
mod camera;
mod fp;
mod hittable;
mod materials;
mod scene;
mod sphere;
mod types;

use acceleration::*;
use camera::*;
use hittable::*;
use rayon::prelude::*;
use rng::*;
use scene::*;
use std::ptr::NonNull;
use std::sync::atomic::{AtomicUsize, Ordering};
use threadpool::*;
use types::*;
use ultraviolet::interp::Lerp;

#[allow(dead_code)]
fn closest_hit<'scene>(
    objects: &[&'scene (dyn hittable::Hittable + Sync + 'scene)],
    r: &Ray,
    t_min: f32,
    t_max: f32,
) -> Option<Hit<'scene>> {
    let mut hit: Option<Hit<'scene>> = None;
    let mut t_max = t_max;
    for o in objects {
        if let Some(h) = o.hit(r, t_min, t_max) {
            t_max = h.t;
            hit = Some(h);
        }
    }
    hit
}

// Traces a ray. This is `color` in the book.
fn trace(
    bvh: &Bvh,
    objects: &[&(dyn Hittable + Sync)],
    r: &Ray,
    depth: i32,
    rng: &mut RttRng,
) -> Vec4f {
    if let Some(hit) = bvh.hit(objects, r, 0.001, std::f32::MAX) {
        if depth >= 50 {
            return Vec4f::zero();
        }
        if let Some(sc) = hit.material.scatter(r, &hit, rng) {
            return sc.attenuation * trace(bvh, objects, &sc.r, depth + 1, rng);
        }
        return Vec4f::zero();
    }
    let white: Vec4f = Vec4f::one();
    let sky_blue: Vec4f = Vec4f::new(0.5, 0.7, 1., 0.);
    let unit_direction = r.direction.normalized();
    let t = 0.5 * (unit_direction.y + 1.);
    white.lerp(sky_blue, t)
}

// Traces a ray bundle.
fn tracex8(
    bvh: &Bvh,
    objects: &[&(dyn Hittable + Sync)],
    r: &Rayx8,
    depth: i32,
    rng: &mut RttRng,
) -> Vec4f {
    if let Some(hit) = bvh.hit(objects, r, 0.001, std::f32::MAX) {
        if depth >= 50 {
            return Vec4f::zero();
        }
        if let Some(sc) = hit.material.scatter(r, &hit, rng) {
            return sc.attenuation * trace(bvh, objects, &sc.r, depth + 1, rng);
        }
        return Vec4f::zero();
    }
    let white: Vec4f = Vec4f::one();
    let sky_blue: Vec4f = Vec4f::new(0.5, 0.7, 1., 0.);
    let unit_direction = r.direction.normalized();
    let t = 0.5 * (unit_direction.y + 1.);
    white.lerp(sky_blue, t)
}

pub fn tracescene(
    nx: usize,
    ny: usize,
    ns: usize,
    scene: &Scene,
    camera: &Camera,
    pool: &rayon::ThreadPool,
    pdc: &AtomicUsize,
) -> Vec<u8> {
    const BYTES_PER_PIXEL: usize = 3;
    let mut pixels = vec![0u8; ny * nx * BYTES_PER_PIXEL];
    let objects = scene.objects();
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
                let mut c = Vec4f::zero();
                for _ in 0..ns {
                    let ray = camera.gen_ray(
                        (x + rng.gen::<f32>()) / nx as f32,
                        (y + rng.gen::<f32>()) / ny as f32,
                        &mut rng,
                    );
                    c += trace(&scene.bvh, &objects, &ray, 0, &mut rng);
                }
                // The book uses a simple gamma 2.0 function, not the sRGB OETF.
                c.apply(|e| (e / (ns as f32)).sqrt() * 255.99);
                chunk[0] = c.x as u8;
                chunk[1] = c.y as u8;
                chunk[2] = c.z as u8;
                pdc.fetch_add(1, Ordering::Relaxed);
            });
    });
    pixels
}
