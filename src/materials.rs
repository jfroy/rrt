use super::hittable::*;
use super::rng::*;
use super::types::*;

pub trait Material {
    fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<ScatteredRay>;
}

// Lambertian

pub struct Lambertian {
    pub albedo: Vec4f,
}

fn near_zero(v: Vec4f) -> bool {
    v.x.abs() < 1e-8 && v.y.abs() < 1e-8 && v.z.abs() < 1e-8 && v.w.abs() < 1e-8
}

impl Material for Lambertian {
    fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<ScatteredRay> {
        let direction = hit.normal + Vec4f::gen_uniform_random_unit(rng);
        Some(ScatteredRay {
            r: Ray {
                origin: hit.p,
                direction: if near_zero(direction) {
                    hit.normal
                } else {
                    direction
                },
                time: r.time,
            },
            attenuation: self.albedo,
        })
    }
}

// Metal

pub struct Metal {
    pub albedo: Vec4f,
    pub fuzz: f32,
}

impl Material for Metal {
    fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<ScatteredRay> {
        let direction = r.direction.normalized().reflected(hit.normal)
            + self.fuzz * Vec4f::gen_uniform_random_in_unit_sphere(rng);
        if direction.dot(hit.normal) > 0. {
            Some(ScatteredRay {
                r: Ray {
                    origin: hit.p,
                    direction,
                    time: r.time,
                },
                attenuation: self.albedo,
            })
        } else {
            None
        }
    }
}

// Dielectric

fn reflectance(cosine: f32, ref_idx: f32) -> f32 {
    // Use Schlick's approximation for reflectance.
    let r0 = ((1. - ref_idx) / (1. + ref_idx)).powf(2.);
    r0 + (1. - r0) * (1. - cosine).powf(5.)
}

pub struct Dielectric {
    pub ref_idx: f32,
}

impl Material for Dielectric {
    fn scatter(&self, r: &Ray, hit: &Hit, rng: &mut RttRng) -> Option<ScatteredRay> {
        let refraction_ratio = if hit.front_face {
            1. / self.ref_idx
        } else {
            self.ref_idx
        };

        let unit_direction = r.direction.normalized();
        let cos_theta = (-unit_direction).dot(hit.normal).min(1.);
        let sin_theta = (1. - (cos_theta * cos_theta)).sqrt();

        let cannot_refract = refraction_ratio * sin_theta > 1.;
        let direction =
            if cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.gen::<f32>() {
                unit_direction.reflected(hit.normal)
            } else {
                unit_direction.refracted(hit.normal, refraction_ratio)
            };

        Some(ScatteredRay {
            r: Ray {
                origin: hit.p,
                direction,
                time: r.time,
            },
            attenuation: Vec4f::new(1., 1., 1., 0.),
        })
    }
}
