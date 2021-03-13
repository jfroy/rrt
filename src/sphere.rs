#![allow(clippy::suspicious_operation_groupings)]

use super::acceleration::*;
use super::materials::*;
use super::types::*;
use std::borrow::Borrow;

pub struct Sphere {
    pub center0: Vec4f,
    pub center1: Vec4f,
    pub time0: f32,
    pub time1: f32,
    pub radius: f32,
    pub material: Box<dyn Material + Sync>,
}

pub struct StationarySphere {
    pub center: Vec4f,
    pub radius: f32,
    pub material: Box<dyn Material + Sync>,
}

impl Sphere {
    pub fn from(s: StationarySphere) -> Sphere {
        Sphere {
            center0: s.center,
            center1: s.center,
            time0: 0.,
            time1: 1.,
            radius: s.radius,
            material: s.material,
        }
    }

    fn center(&self, t: f32) -> Vec4f {
        self.center0
            + ((t - self.time0) / (self.time1 - self.time0)) * (self.center1 - self.center0)
    }
}

impl Hittable for Sphere {
    fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>> {
        let center = self.center(r.time);
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
        let p = r.point_at_parameter(root);
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

    fn aabb(&self, time0: f32, time1: f32) -> Option<Aabb> {
        let box0 = Aabb {
            minimum: self.center(time0) - Vec4f::broadcast(self.radius).with_w(0.),
            maximum: self.center(time0) + Vec4f::broadcast(self.radius).with_w(0.),
        };
        let box1 = Aabb {
            minimum: self.center(time1) - Vec4f::broadcast(self.radius).with_w(0.),
            maximum: self.center(time1) + Vec4f::broadcast(self.radius).with_w(0.),
        };
        Some(Aabb::surrounding(box0, box1))
    }
}
