#![allow(clippy::suspicious_operation_groupings)]

use super::acceleration::*;
use super::hittable::*;
use super::materials::*;
use super::types::*;
use std::borrow::Borrow;
use ultraviolet as uv;
use wide::*;

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
        // NOTE: rrt does not normalize ray directions.
        let a = r.direction.mag_sq();
        let half_b = oc.dot(r.direction);
        let c = oc.mag_sq() - (self.radius * self.radius);
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

struct Spheresx8 {
    center0: uv::Vec4x8,
    center1: uv::Vec4x8,
    time0: uv::f32x8,
    time1: uv::f32x8,
    radius: uv::f32x8,
}

impl Spheresx8 {
    fn center_at(&self, t: uv::f32x8) -> uv::Vec4x8 {
        self.center0
            + ((t - self.time0) / (self.time1 - self.time0)) * (self.center1 - self.center0)
    }
}

macro_rules! pack_spheres_x8 {
    ( $o:ident, $v:ident, $( $c:ident ),+ ) => {{
        (
            $({
                let mut f = [0.; 8];
                for i in 0..8 {
                    f[i] = $o[i].$v.$c;
                };
                f.into()
            }),*
        ).into()}
    };
    ( $o:ident, $v:ident ) => {{
        let mut f = [0.; 8];
        for i in 0..8 {
            f[i] = $o[i].$v;
        };
        f.into()
    }};
}

impl From<[&Sphere; 8]> for Spheresx8 {
    fn from(spheres: [&Sphere; 8]) -> Self {
        let center0: uv::Vec4x8 = pack_spheres_x8!(spheres, center0, x, y, z, w);
        let center1: uv::Vec4x8 = pack_spheres_x8!(spheres, center1, x, y, z, w);
        let time0: uv::f32x8 = pack_spheres_x8!(spheres, time0);
        let time1: uv::f32x8 = pack_spheres_x8!(spheres, time1);
        let radius: uv::f32x8 = pack_spheres_x8!(spheres, radius);

        Spheresx8 {
            center0,
            center1,
            time0,
            time1,
            radius,
        }
    }
}

pub fn ray_sphere_intersect_x8<'scene>(
    spheres: [&'scene Sphere; 8],
    r: &Rayx8,
    t_min: uv::f32x8,
    t_max: uv::f32x8,
) -> Hitx8<'scene> {
    let ps: Spheresx8 = spheres.into();
    let center = ps.center_at(r.time);
    let oc = r.origin - center;
    // NOTE: rrt does not normalize ray directions.
    let a = r.direction.mag_sq();
    let half_b = oc.dot(r.direction);
    let c = oc.mag_sq() - (ps.radius * ps.radius);
    let discriminant = (half_b * half_b) - (a * c);
    let discriminant_not_neg_mask = discriminant.cmp_ge(wide::f32x8::ZERO);
    let discriminant_sqrt = discriminant.sqrt();
    let root1 = (-half_b - discriminant_sqrt) / a;
    let root2 = (-half_b + discriminant_sqrt) / a;
    let root1_valid = root1.cmp_ge(t_min) & root1.cmp_le(t_max) & discriminant_not_neg_mask;
    let root2_valid = root2.cmp_ge(t_min) & root2.cmp_le(t_max) & discriminant_not_neg_mask;
    let root = root1_valid.blend(root1, root2);
    let p = r.point_at(root);
    let outward_normal = (p - center) / ps.radius;
    let front_face_mask = r.direction.dot(outward_normal).cmp_lt(wide::f32x8::ZERO);
    let normal = uv::Vec4x8::blend(front_face_mask, outward_normal, -outward_normal);
    Hitx8 {
        p,
        normal,
        t: root,
        material: [
            Some(spheres[0].material.borrow()),
            Some(spheres[1].material.borrow()),
            Some(spheres[2].material.borrow()),
            Some(spheres[3].material.borrow()),
            Some(spheres[4].material.borrow()),
            Some(spheres[5].material.borrow()),
            Some(spheres[6].material.borrow()),
            Some(spheres[7].material.borrow()),
        ],
        front_face_bmask: front_face_mask.move_mask(),
        valid_bmask: (root1_valid | root2_valid).move_mask(),
    }
}
