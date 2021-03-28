use super::acceleration::*;
use super::materials::*;
use super::types::*;

pub trait Hittable {
    fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>>;
    fn aabb(&self) -> Aabb;
}

pub struct Hit<'scene> {
    pub p: Vec4f,
    pub normal: Vec4f,
    pub t: f32,
    pub material: &'scene (dyn Material + Sync),
    pub front_face: bool,
}

#[derive(Clone, Copy)]
pub struct Unhittable {}

impl Hittable for Unhittable {
    fn hit<'scene>(&'scene self, _: &Ray, _: f32, _: f32) -> Option<Hit<'scene>> {
        None
    }
    fn aabb(&self) -> Aabb {
        Aabb::default()
    }
}
