use rand::SeedableRng;
use rrt::{chap12_scene, init_pool_with_rng, tracescene, RttRng};

fn main() {
  const NX: usize = 1200;
  const NY: usize = 800;
  const NS: usize = 10;

  eprintln!("Rendering {} x {} image using {} samples.", NX, NY, NS);

  let mut rng = RttRng::seed_from_u64(0);
  let (scene, camera) = chap12_scene(NX, NY, &mut rng);
  let pool = init_pool_with_rng(rng);
  let pixels = tracescene(NX, NY, NS, &scene, &camera, &pool);
  image::save_buffer(
    "o.ppm",
    &pixels[..],
    NX as u32,
    NY as u32,
    image::ColorType::Rgb8,
  )
  .unwrap()
}
