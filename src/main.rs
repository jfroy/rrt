#![warn(clippy::all)]

use rrt::tracescene;

fn main() {
  const NX: usize = 1000;
  const NY: usize = 500;
  const NS: usize = 100;
  eprintln!("Rendering {} x {} image using {} samples.", NX, NY, NS);
  let pixels = tracescene(NX, NY, NS);
  image::save_buffer("o.ppm", &pixels[..], NX as u32, NY as u32, image::RGB(8)).unwrap()
}
