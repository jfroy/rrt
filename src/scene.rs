use super::acceleration::*;
use super::hittable::*;
use super::rng::*;
use super::sphere::*;
use super::types::*;
use rand::distributions::Standard;

pub struct Scene<'a> {
    pub spheres: Vec<Sphere>,
    pub bvh_nodes: Vec<BvhNode<'a>>,
}

impl<'a> Scene<'a> {
    pub fn new() -> Scene<'a> {
        Scene {
            spheres: vec![],
            bvh_nodes: vec![],
        }
    }

    pub fn build_bvh(&self, rng: &mut RttRng) {
        let mut nodes: Vec<BvhNode> = Vec::new();
        let mut refs: Vec<&Sphere> = self.spheres.iter().collect();
        let mut slices: Vec<&mut [&Sphere]> = Vec::new();
        slices.push(&mut refs[..]);
        while let Some(s) = slices.pop() {
            let axis: Axis = rng.sample(Standard);
            match s.len() {
                0 => {}
                1 => nodes.push(BvhNode::new(s[0], s[0]).unwrap()),
                2 => nodes.push(BvhNode::new(s[0], s[1]).unwrap()),
                _ => {
                    s.sort_unstable_by(|a, b| a.aabb.axis_cmp(&b.aabb, axis));
                    let (left, right) = s.split_at_mut(s.len() / 2);
                    slices.push(left);
                    slices.push(right);
                }
            }
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
