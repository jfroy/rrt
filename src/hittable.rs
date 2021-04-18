use super::acceleration::*;
use super::materials::*;
use super::types::*;
use ultraviolet as uv;

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
#[derive(Default)]
pub struct Hitx8<'scene> {
    pub p: uv::Vec4x8,
    pub normal: uv::Vec4x8,
    pub t: uv::f32x8,
    pub material: [Option<&'scene (dyn Material + Sync)>; 8],
    pub front_face_bmask: i32,
    pub valid_bmask: i32,
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
