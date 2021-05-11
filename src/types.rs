// Provides basic types.
#![allow(dead_code)]

use ultraviolet as uv;

pub type Vec3f = uv::Vec3;
pub type Vec4f = uv::Vec4;

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
pub struct Rayx8 {
    pub origin: uv::Vec4x8,
    pub direction: uv::Vec4x8,
    pub time: uv::f32x8,
}

impl Rayx8 {
    pub fn point_at(&self, t: uv::f32x8) -> uv::Vec4x8 {
        self.origin + (t * self.direction)
    }

    pub fn inv_direction(&self) -> uv::Vec4x8 {
        uv::Vec4x8::one() / self.direction
    }
}

#[derive(Clone, Copy)]
pub struct ScatteredRayx8 {
    pub r: Rayx8,
    pub attenuation: uv::Vec4x8,
}

#[derive(Clone, Copy)]
pub enum Axis {
    X = 0,
    Y,
    Z,
}

pub type Axisx8 = [Axis; 8];
