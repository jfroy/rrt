#![allow(clippy::suspicious_operation_groupings)]

use super::acceleration::*;
use super::hittable::*;
use super::materials::*;
use super::types::*;
use std::borrow::Borrow;

pub struct Sphere {
    center0: Vec4f,
    center1: Vec4f,
    time0: f32,
    time1: f32,
    radius: f32,
    material: Box<dyn Material + Sync>,
    aabb: Aabb,
}

pub struct StationarySphere {
    pub center: Vec4f,
    pub radius: f32,
    pub material: Box<dyn Material + Sync>,
}

fn sphere_aabb(center0: Vec4f, center1: Vec4f, radius: f32) -> Aabb {
    let box0 = Aabb {
        minimum: center0 - Vec4f::broadcast(radius),
        maximum: center0 + Vec4f::broadcast(radius),
    };
    let box1 = Aabb {
        minimum: center1 - Vec4f::broadcast(radius),
        maximum: center1 + Vec4f::broadcast(radius),
    };
    box0.union(box1)
}

impl Sphere {
    pub fn new(
        center0: Vec4f,
        center1: Vec4f,
        time0: f32,
        time1: f32,
        radius: f32,
        material: Box<dyn Material + Sync>,
    ) -> Sphere {
        Sphere {
            center0,
            center1,
            time0,
            time1,
            radius,
            material,
            aabb: sphere_aabb(center0, center1, radius),
        }
    }

    pub fn from(s: StationarySphere) -> Sphere {
        Sphere::new(s.center, s.center, 0., 1., s.radius, s.material)
    }

    fn center_at(&self, t: f32) -> Vec4f {
        self.center0
            + ((t - self.time0) / (self.time1 - self.time0)) * (self.center1 - self.center0)
    }
}

impl Hittable for Sphere {
    fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>> {
        let center = self.center_at(r.time);
        let oc = r.origin - center;
        let a = r.direction.dot(r.direction);
        let half_b = oc.dot(r.direction);
        let c = oc.dot(oc) - (self.radius * self.radius);
        let discriminant = (half_b * half_b) - (a * c);
        if discriminant < 0. {
            return None;
        }
        let discriminant_sqrt = discriminant.sqrt();
        let mut root = (-half_b - discriminant_sqrt) / a;
        if root < t_min || t_max < root {
            root = (-half_b + discriminant_sqrt) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }
        let p = r.point_at(root);
        let outward_normal = (p - center) / self.radius;
        let front_face = r.direction.dot(outward_normal) < 0.;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        Some(Hit {
            p,
            normal,
            t: root,
            material: self.material.borrow(),
            front_face,
        })
    }

    fn aabb(&self) -> Aabb {
        self.aabb
    }
}
