// Provides basic types.
#![allow(dead_code)]

use ultraviolet;

pub type Vec3f = ultraviolet::Vec3;
pub type Vec4f = ultraviolet::Vec4;

#[derive(Clone, Copy)]
pub struct Ray {
    pub origin: Vec4f,
    pub direction: Vec4f,
    pub time: f32,
}

impl Ray {
    pub fn point_at(&self, t: f32) -> Vec4f {
        self.origin + (t * self.direction)
    }

    pub fn inv_direction(&self) -> Vec4f {
        Vec4f::one() / self.direction
    }
}

#[derive(Clone, Copy)]
pub struct ScatteredRay {
    pub r: Ray,
    pub attenuation: Vec4f,
}

#[derive(Clone, Copy)]
pub enum Axis {
    X = 0,
    Y,
    Z,
}
