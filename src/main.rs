use clap::{App, Arg};
use pbr::ProgressBar;
use rrt::book2chap2::*;
use rrt::rng::*;
use rrt::threadpool::*;
use rrt::tracescene;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{thread, time};

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
        .arg(
            Arg::new("seed")
                .long("seed")
                .short('e')
                .takes_value(true)
                .default_value("0")
                .about("rng seed"),
        )
        .arg(
            Arg::new("random")
                .long("random")
                .short('m')
                .about("use a random rng seed"),
        )
        .get_matches();
    let (w, h) = match parse_resolution(arg_matches.value_of("resolution")) {
        Some(v) => v,
        None => {
            eprintln!("invalid resolution");
            return;
        }
    };
    let spp: usize = match parse_arg(arg_matches.value_of("samples")) {
        Some(v) => v,
        None => {
            eprintln!("invalid sample count");
            return;
        }
    };
    let mut rng: RttRng = if arg_matches.is_present("random") {
        RttRng::from_entropy()
    } else {
        match parse_arg(arg_matches.value_of("seed")) {
            Some(v) => RttRng::seed_from_u64(v),
            None => {
                eprintln!("invalid rng seed");
                return;
            }
        }
    };

    eprintln!(
        "Rendering {} x {} image using {} samples per pixel.",
        w, h, spp
    );

    let pxcount = Arc::new(AtomicUsize::new(0));
    let ui_pxcount = Arc::clone(&pxcount);
    let ui_thread = thread::Builder::new()
        .name("ui".to_string())
        .spawn(move || {
            let t = w * h;
            let mut pb = ProgressBar::new((t) as u64);
            loop {
                let x = ui_pxcount.load(Ordering::Relaxed);
                pb.set(x as u64);
                thread::sleep(time::Duration::from_secs(1));
                if x >= t {
                    break;
                }
            }
        })
        .unwrap();

    let (scene, camera) = book2_chap2_scene(w, h, &mut rng);
    let pool = init_pool_with_rng(rng);
    let pixels = tracescene(w, h, spp, &scene, &camera, &pool, &pxcount);
    ui_thread.join().unwrap();
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

fn parse_arg<T: std::str::FromStr>(s: Option<&str>) -> Option<T> {
    match s?.parse::<T>() {
        Ok(n) => Some(n),
        Err(_) => None,
    }
}
