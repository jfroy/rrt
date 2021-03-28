use crate::fp::gamma_eb;
use crate::hittable::*;
use crate::simd_llvm;
use crate::types::*;
use arrayvec::ArrayVec;
use itertools;
use std::cmp::Ordering;

#[derive(Clone, Copy)]
pub struct Aabb {
    pub minimum: Vec4f,
    pub maximum: Vec4f,
}

impl Aabb {
    pub fn union(self: Aabb, b: Aabb) -> Aabb {
        let minimum = Vec4f::partial_min(self.minimum, b.minimum);
        let maximum = Vec4f::partial_max(self.maximum, b.maximum);
        Aabb { minimum, maximum }
    }

    pub fn union_p(self: Aabb, p: Vec4f) -> Aabb {
        let minimum = Vec4f::partial_min(self.minimum, p);
        let maximum = Vec4f::partial_max(self.maximum, p);
        Aabb { minimum, maximum }
    }

    pub fn centroid(self) -> Vec4f {
        0.5 * self.minimum + 0.5 * self.maximum
    }

    pub fn diagonal(self) -> Vec4f {
        self.maximum - self.minimum
    }

    pub fn max_extent(self) -> Axis {
        let d = self.diagonal();
        if d.x > d.y && d.x > d.z {
            Axis::X
        } else if d.y > d.z {
            Axis::Y
        } else {
            Axis::Z
        }
    }

    pub fn axis_cmp(self, other: Aabb, axis: Axis) -> Ordering {
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

impl Default for Aabb {
    fn default() -> Self {
        Aabb {
            minimum: Vec4f::broadcast(f32::MAX),
            maximum: Vec4f::broadcast(f32::MIN),
        }
    }
}

pub struct Bvh {
    nodes: Vec<BvhNode>,
}

#[derive(Clone, Copy)]
pub enum BvhPartitionMethod {
    Middle,
    EqualCount,
}

enum BvhNode {
    Leaf {
        obj_i: usize,
    },
    Inner {
        aabb: Aabb,
        axis: Axis,
        left_i: Option<usize>,
        right_i: Option<usize>,
    },
}

struct BvhObjectInfo<'a> {
    pub index: usize,
    pub object: &'a (dyn Hittable + Sync),
    pub aabb: Aabb,
    pub centroid: Vec4f,
}

impl Bvh {
    pub fn new() -> Bvh {
        Bvh { nodes: vec![] }
    }

    /// Build the BVH using a slice of objects and a partitioning method.
    pub fn build(&mut self, objects: &[&(dyn Hittable + Sync)], method: BvhPartitionMethod) {
        let mut infos: Vec<_> = objects
            .iter()
            .enumerate()
            .map(|(index, object)| {
                let aabb = object.aabb();
                BvhObjectInfo {
                    index,
                    object: *object,
                    aabb,
                    centroid: aabb.centroid(),
                }
            })
            .collect();
        self.build_r(&mut infos, method);
    }

    /// Build the BVH recursively using a slice of objects and a partitioning
    /// method. See http://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies.html#BVHConstruction
    fn build_r(
        &mut self,
        infos: &mut [BvhObjectInfo],
        mut method: BvhPartitionMethod,
    ) -> Option<usize> {
        match infos.len() {
            0 => None,
            1 => {
                self.nodes.push(BvhNode::Leaf {
                    obj_i: infos[0].index,
                });
                Some(self.nodes.len() - 1)
            }
            _ => {
                let aabb = infos
                    .iter()
                    .fold(Aabb::default(), |aabb, o| aabb.union(o.aabb));
                let centroid_bounds = infos
                    .iter()
                    .fold(Aabb::default(), |aabb, o| aabb.union_p(o.centroid));
                let axis = centroid_bounds.max_extent();
                let (left_o, right_o) = loop {
                    let (left_o, right_o) = match method {
                        BvhPartitionMethod::Middle => {
                            let midpoint = (centroid_bounds.minimum[axis as usize]
                                + centroid_bounds.maximum[axis as usize])
                                / 2.;
                            let pivot_i = itertools::partition(infos.iter_mut(), |o| {
                                o.centroid[axis as usize] < midpoint
                            });
                            if pivot_i == 0 || pivot_i == infos.len() {
                                method = BvhPartitionMethod::EqualCount;
                                continue;
                            }
                            infos.split_at_mut(pivot_i)
                        }
                        BvhPartitionMethod::EqualCount => {
                            let pivot_i = infos.len() / 2;
                            infos.select_nth_unstable_by(pivot_i, |a, b| {
                                a.centroid[axis as usize]
                                    .partial_cmp(&b.centroid[axis as usize])
                                    .unwrap_or(Ordering::Equal)
                            });
                            infos.split_at_mut(pivot_i)
                        }
                    };
                    break (left_o, right_o);
                };
                self.nodes.push(BvhNode::Inner {
                    aabb,
                    axis,
                    left_i: None,
                    right_i: None,
                });
                let inner_i = self.nodes.len() - 1;
                let new_left_i = self.build_r(left_o, method);
                let new_right_i = self.build_r(right_o, method);
                if let BvhNode::Inner {
                    left_i, right_i, ..
                } = &mut self.nodes[inner_i]
                {
                    *left_i = new_left_i;
                    *right_i = new_right_i;
                }
                Some(inner_i)
            }
        }
    }

    /// Hit traverses the BVH and returns the closest hit, if any.
    /// See http://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies.html#Traversal
    pub fn hit<'scene>(
        &self,
        objects: &[&'scene (dyn Hittable + Sync)],
        r: &Ray,
        t_min: f32,
        t_max: f32,
    ) -> Option<Hit<'scene>> {
        let inv_dir = r.direction.with_w(1.).recip();
        let inv_dir_neg: Vec4i8 = inv_dir.partial_cmplt(&Vec4f::zero()).as_();
        let mut index_stack = ArrayVec::<_, 64>::new();
        index_stack.push(0);
        let mut hit: Option<Hit<'scene>> = None;
        let mut t_max = t_max;
        while let Some(i) = index_stack.pop() {
            match &self.nodes[i] {
                BvhNode::Leaf { obj_i } => {
                    if let Some(h) = objects[*obj_i].hit(r, t_min, t_max) {
                        t_max = h.t;
                        hit = Some(h);
                    }
                }
                BvhNode::Inner {
                    aabb,
                    axis,
                    left_i,
                    right_i,
                } => {
                    if !aabb.intersect(r, inv_dir, inv_dir_neg, t_min, t_max) {
                        continue;
                    }
                    if inv_dir_neg[*axis as usize] != 0 {
                        if let Some(i) = *left_i {
                            index_stack.push(i);
                        }
                        if let Some(i) = *right_i {
                            index_stack.push(i);
                        }
                    } else {
                        if let Some(i) = *right_i {
                            index_stack.push(i);
                        }
                        if let Some(i) = *left_i {
                            index_stack.push(i);
                        }
                    }
                }
            }
        }
        hit
    }
}
