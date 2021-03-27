use crate::fp::gamma_eb;
use crate::hittable::*;
use crate::simd_llvm;
use crate::types::*;
use arrayvec::ArrayVec;
use std::cmp::Ordering;

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

    pub fn axis_cmp(&self, other: &Aabb, axis: Axis) -> Ordering {
        let a = self.minimum[axis as usize];
        let b = other.minimum[axis as usize];
        a.partial_cmp(&b).unwrap_or(Ordering::Equal)
    }

    /// Intersect returns true if a ray conservatively intersects the AABB.
    ///
    /// This implementation is essentially copied from [pbrt]. See also the
    /// section on [conservative intersections].
    ///
    /// [pbrt]: http://www.pbr-book.org/3ed-2018/Shapes/Basic_Shape_Interface.html#RayndashBoundsIntersections
    /// [conservative intersections]: http://www.pbr-book.org/3ed-2018/Shapes/Managing_Rounding_Error.html#ConservativeRayndashBoundsIntersections
    pub fn intersect(&self, r: &Ray, inv_d: Vec4f, inv_dir_neg: Vec4i8, t0: f32, t1: f32) -> bool {
        let min: Vec4f;
        let max: Vec4f;
        unsafe {
            min = simd_llvm::simd_select(inv_dir_neg, self.maximum, self.minimum);
            max = simd_llvm::simd_select(inv_dir_neg, self.minimum, self.maximum);
        }
        let t_near = (min - r.origin) * inv_d;
        let t_far = (max - r.origin) * inv_d * (1. + 2. * gamma_eb(3));
        let t0 = Vec4f::reduce_partial_max(t_near.with_w(t0));
        let t1 = Vec4f::reduce_partial_min(t_far.with_w(t1));
        t0 <= t1
    }
}

pub struct Bvh<'a> {
    nodes: Vec<BvhNode<'a>>,
}

impl<'a> Bvh<'a> {
    pub fn new() -> Bvh<'a> {
        Bvh { nodes: vec![] }
    }
}

pub enum BvhNode<'a> {
    Leaf {
        hittable: &'a (dyn Hittable + Sync),
    },
    Inner {
        aabb: Aabb,
        axis: Axis,
        left: usize,
        right: usize,
    },
}

impl<'a> Hittable for Bvh<'a> {
    /// Hit traverses the BVH and returns the closest hit, if any.
    /// See http://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies.html#Traversal
    fn hit(&self, r: &Ray, t_min: f32, t_max: f32) -> Option<Hit<'a>> {
        let inv_dir = r.direction.with_w(1.).recip();
        let inv_dir_neg: Vec4i8 = inv_dir.partial_cmplt(&Vec4f::zero()).as_();
        let mut node_index_stack = ArrayVec::<_, 64>::new();
        let mut hit: Option<Hit<'a>> = None;
        let mut t_max = t_max;
        while let Some(i) = node_index_stack.pop() {
            let node = &self.nodes[i];
            match node {
                BvhNode::Leaf { hittable } => {
                    if let Some(h) = hittable.hit(r, t_min, t_max) {
                        t_max = h.t;
                        hit = Some(h);
                    }
                }
                BvhNode::Inner {
                    aabb,
                    axis,
                    left,
                    right,
                } => {
                    if !aabb.intersect(r, inv_dir, inv_dir_neg, t_min, t_max) {
                        continue;
                    }
                    if inv_dir_neg[*axis as usize] != 0 {
                        node_index_stack.push(*left);
                        node_index_stack.push(*right);
                    } else {
                        node_index_stack.push(*right);
                        node_index_stack.push(*left);
                    }
                }
            }
        }
        hit
    }

    fn aabb(&self) -> Option<Aabb> {
        None
    }
}
