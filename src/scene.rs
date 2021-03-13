use super::acceleration::*;
use super::materials::*;
use super::sphere::*;
use super::types::*;

pub struct Scene {
    pub spheres: Vec<Sphere>,
}

impl Scene {
    pub fn new() -> Scene {
        Scene { spheres: vec![] }
    }
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

    fn aabb(&self, time0: f32, time1: f32) -> Option<Aabb> {
        if self.spheres.is_empty() {
            return None;
        }
        let mut bb = Aabb::zero();
        for sphere in &self.spheres {
            if let Some(aabb) = sphere.aabb(time0, time1) {
                bb = Aabb::surrounding(bb, aabb);
            } else {
                return None;
            }
        }
        Some(bb)
    }
}
