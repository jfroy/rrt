use clap::{App, Arg};
use rrt::chap12::*;
use rrt::rng::*;
use rrt::threadpool::*;
use rrt::tracescene;

fn main() {
  let arg_matches = App::new("rrt")
    .version("0.1.0")
    .author("Jean-Francois Roy <jf@devklog.net>")
    .about("Ray Tracing In One Weekend")
    .arg(
      Arg::new("resolution")
        .long("resolution")
        .short('r')
        .takes_value(true)
        .default_value("1200x800")
        .about("output resolution in pixels"),
    )
    .arg(
      Arg::new("samples")
        .long("samples")
        .short('s')
        .takes_value(true)
        .default_value("10")
        .about("samples per pixel"),
    )
    .get_matches();
  let (w, h) = match parse_resolution(arg_matches.value_of("resolution")) {
    Some(v) => v,
    None => {
      eprintln!("invalid resolution");
      return;
    }
  };
  let samples = match parse_samples(arg_matches.value_of("samples")) {
    Some(v) => v,
    None => {
      eprintln!("invalid sample count");
      return;
    }
  };

  eprintln!("Rendering {} x {} image using {} samples.", w, h, samples);

  let mut rng = RttRng::seed_from_u64(0);
  let (scene, camera) = chap12_scene(w, h, &mut rng);
  let pool = init_pool_with_rng(rng);
  let pixels = tracescene(w, h, samples, &scene, &camera, &pool);
  image::save_buffer(
    "o.ppm",
    &pixels[..],
    w as u32,
    h as u32,
    image::ColorType::Rgb8,
  )
  .unwrap()
}

fn parse_resolution(s: Option<&str>) -> Option<(usize, usize)> {
  let v: Vec<&str> = s?.split('x').collect();
  if v.len() != 2 {
    return None;
  }
  let w = match v[0].parse::<usize>() {
    Ok(n) => n,
    Err(_) => {
      return None;
    }
  };
  let h = match v[1].parse::<usize>() {
    Ok(n) => n,
    Err(_) => {
      return None;
    }
  };
  Some((w, h))
}

fn parse_samples(s: Option<&str>) -> Option<usize> {
  match s?.parse::<usize>() {
    Ok(n) => Some(n),
    Err(_) => None,
  }
}
