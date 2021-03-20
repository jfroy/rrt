use super::acceleration::*;
use super::materials::*;
use super::sphere::*;
use super::types::*;

pub struct Scene<'a> {
    pub spheres: Vec<Sphere>,
    pub bvh_nodes: Vec<Bvh<'a>>,
}

impl<'a> Scene<'a> {
    pub fn new() -> Scene<'a> {
        Scene {
            spheres: vec![],
            bvh_nodes: vec![],
        }
    }
}

impl<'a> Hittable for Scene<'a> {
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

    fn aabb(&self) -> Option<Aabb> {
        if self.spheres.is_empty() {
            return None;
        }
        let mut bb = Aabb::zero();
        for sphere in &self.spheres {
            if let Some(aabb) = sphere.aabb() {
                bb = Aabb::surrounding(bb, aabb);
            } else {
                return None;
            }
        }
        Some(bb)
    }
}
