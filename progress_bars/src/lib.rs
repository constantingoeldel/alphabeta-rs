use std::ops::Deref;

use indicatif::ProgressStyle;
pub use indicatif::{MultiProgress, ProgressBar};

pub struct Progress(pub ProgressBar);

pub fn progress_style() -> ProgressStyle {
    ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
        .unwrap()
        .progress_chars("##-")
}

impl Progress {
    pub fn new(name: &'static str, iterations: usize) -> Self {
        let pb = ProgressBar::new(iterations as u64);

        pb.set_message(name);
        pb.set_style(
            ProgressStyle::with_template(
                "{msg} [{elapsed}] {wide_bar:40.cyan/blue} {pos:>7}/{len:7}",
            )
            .unwrap(),
        );
        Progress(pb)
    }
}

impl Deref for Progress {
    type Target = ProgressBar;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
