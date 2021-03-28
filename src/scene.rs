use super::acceleration::*;
use super::hittable::*;
use super::sphere::*;

pub struct Scene {
    pub spheres: Vec<Sphere>,
    pub bvh: Bvh,
}

fn as_object_refs<H>(objects: &[H]) -> Vec<&(dyn Hittable + Sync)>
where
    H: Hittable + Sync,
{
    objects
        .iter()
        .map(|e| e as &(dyn Hittable + Sync))
        .collect()
}

impl Scene {
    pub fn new() -> Scene {
        Scene {
            spheres: vec![],
            bvh: Bvh::new(),
        }
    }

    pub fn objects(&self) -> Vec<&(dyn Hittable + Sync)> {
        as_object_refs(&self.spheres)
    }

    pub fn build_bvh(&mut self, method: BvhPartitionMethod) {
        let objects = as_object_refs(&self.spheres);
        self.bvh.build(&objects[..], method);
    }
}
