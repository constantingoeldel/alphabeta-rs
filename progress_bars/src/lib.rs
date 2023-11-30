use std::{fmt::Write, time::Duration};
use std::ops::Deref;

use indicatif::{HumanDuration, ProgressState, ProgressStyle};
pub use indicatif::{MultiProgress, ProgressBar};

pub struct Progress(pub ProgressBar);

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

pub fn specific(bars: &MultiProgress, iterations: usize) -> (ProgressBar, ProgressBar) {
    let pb_neutral = bars.insert(0, ProgressBar::new(iterations as u64));
    let pb_boot = bars.insert(1, ProgressBar::new(iterations as u64));

    pb_neutral.set_message("ABNeutral");
    pb_boot.set_message("BootModel");
    pb_neutral.set_style(
        ProgressStyle::with_template("{msg} {bar:40.cyan/blue} [{elapsed}] {pos:>7}/{len:7}")
            .unwrap(),
    );

    pb_boot.set_style(
        ProgressStyle::with_template("{msg} {bar:40.cyan/blue} [{elapsed}] {pos:>7}/{len:7}")
            .unwrap(),
    );
    pb_boot.tick();
    (pb_neutral, pb_boot)
}
pub fn multi(total_steps: u64) -> (MultiProgress, ProgressBar) {
    let multi = MultiProgress::new();
    let pb = multi.add(ProgressBar::new(total_steps));
    pb.set_message("Progress ");
    pb.enable_steady_tick(Duration::new(1, 0));
    pb.set_style(
        ProgressStyle::with_template(
            "{msg} {bar:40.magenta/blue} [{elapsed}] {pos:>7}/{len:7} ETA: {eta}",
        )
        .unwrap()
        .with_key("eta", |state: &ProgressState, w: &mut dyn Write| {
            write!(
                w,
                "{}",
                HumanDuration(Duration::from_secs(state.eta().as_secs()))
            )
            .unwrap();
        }),
    );
    (multi, pb)
}
