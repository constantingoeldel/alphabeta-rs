use std::{fmt::Write, time::Duration};

use indicatif::{HumanDuration, MultiProgress, ProgressBar, ProgressState, ProgressStyle};

pub fn specific(bars: &MultiProgress, iterations: u64) -> (ProgressBar, ProgressBar) {
    let pb_neutral = bars.add(ProgressBar::new(iterations));
    let pb_boot = bars.add(ProgressBar::new(iterations));

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
