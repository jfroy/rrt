use super::hittable::*;
use super::types::*;

#[derive(Clone, Copy)]
pub struct Aabb {
    pub minimum: Vec4f,
    pub maximum: Vec4f,
}

impl Aabb {
    pub fn surrounding(a: Aabb, b: Aabb) -> Aabb {
        let minimum = Vec4f::partial_min(a.minimum, b.minimum);
        let maximum = Vec4f::partial_max(a.maximum, b.maximum);
        Aabb { minimum, maximum }
    }

    pub fn zero() -> Aabb {
        Aabb {
            minimum: Vec4f::zero(),
            maximum: Vec4f::zero(),
        }
    }

    pub fn minmax() -> Aabb {
        Aabb {
            minimum: Vec4f::broadcast(f32::MIN).with_w(0.),
            maximum: Vec4f::broadcast(f32::MAX).with_w(0.),
        }
    }

    pub fn hit(&self, r: &Ray, t_min: f32, t_max: f32) -> bool {
        let inv_d = r.direction.with_w(1.).recip();
        let t0 = (self.minimum - r.origin) * inv_d;
        let t1 = (self.maximum - r.origin) * inv_d;
        let t_min = Vec4f::partial_max(t0, Vec4f::broadcast(t_min));
        let t_max = Vec4f::partial_min(t1, Vec4f::broadcast(t_max));
        t_max.partial_cmpgt(&t_min).reduce_and()
    }
}

pub struct Bvh<'a> {
    pub left: &'a (dyn Hittable + Sync),
    pub right: &'a (dyn Hittable + Sync),
    pub aabb: Aabb,
}

impl<'a> Hittable for Bvh<'a> {
    fn hit<'scene>(&'scene self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'scene>> {
        if !self.aabb.hit(r, t_min, t_max) {
            return None;
        }
        let lh = self.left.hit(r, t_min, t_max);
        let rh = self
            .right
            .hit(r, t_min, if let Some(ref hit) = lh { hit.t } else { t_max });
        if rh.is_some() {
            rh
        } else {
            lh
        }
    }

    fn aabb(&self) -> Option<Aabb> {
        Some(self.aabb)
    }
}
