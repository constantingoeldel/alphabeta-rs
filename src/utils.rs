use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

pub fn progress_bar(bars: &MultiProgress, iterations: &u64) -> (ProgressBar, ProgressBar) {
    let pb_neutral = bars.insert(0, ProgressBar::new(*iterations));
    let pb_boot = bars.insert(1, ProgressBar::new(*iterations));

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
