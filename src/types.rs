// Provides basic types.
#![allow(dead_code)]

use vek::vec::repr_simd::*;

pub type Vec3f = vec3::Vec3<f32>;
pub type Vec4b = vec4::Vec4<bool>;
pub type Vec4i8 = vec4::Vec4<i8>;
pub type Vec4f = vec4::Vec4<f32>;
pub type Rgbf32 = rgb::Rgb<f32>;

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
        self.direction.recip().with_w(0.)
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
