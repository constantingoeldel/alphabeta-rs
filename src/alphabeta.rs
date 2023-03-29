use anyhow::anyhow;
use indicatif::MultiProgress;

use crate::{
    arguments,
    pedigree::Pedigree,
    structs::{Model, StandardDeviations},
    utils::progress_bar,
    *,
};

pub struct AlphaBeta<'bar> {
    args: arguments::AlphaBeta,
    bars: &'bar MultiProgress,
}

impl<'bar> AlphaBeta<'bar> {
    pub fn new(args: arguments::AlphaBeta, bars: &'bar MultiProgress) -> Self {
        AlphaBeta { args, bars }
    }

    pub fn run(&self) -> Result<(Model, StandardDeviations)> {
        let arguments::AlphaBeta {
            nodelist,
            edgelist,
            iterations,
            output,
            posterior_max_filter,
        } = &self.args;

        let (pedigree, p0uu) = Pedigree::build(nodelist, edgelist, *posterior_max_filter)
            .map_err(|e| anyhow!("Error while building pedigree: {}", e))?;
        pedigree.to_file(output)?;

        let (pb_neutral, pb_boot) = progress_bar(self.bars, iterations);

        let model = ABneutral::run(
            &pedigree,
            p0uu,
            p0uu,
            1.0,
            *iterations,
            Some(pb_neutral.clone()),
        )
        .map_err(|e| anyhow!("Model failed: {}", e))?;
        let result = BootModel::run(
            &pedigree,
            &model,
            p0uu,
            p0uu,
            1.0,
            *iterations,
            Some(pb_boot.clone()),
        )
        .map_err(|e| anyhow!("Bootstrap failed: {}", e))?;
        self.bars.remove(&pb_neutral);
        self.bars.remove(&pb_boot);

        model.to_file(output, &result)?;
        Ok((model, result))
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{assert_close, assert_close_10_percent, pedigree::Pedigree, ABneutral, BootModel};

    #[test] // Recommended to run with --release
    fn end_to_end() {
        let (pedigree, p0uu) = Pedigree::build(
            Path::new("./data/desired_output/nodelist.txt"),
            Path::new("./data/desired_output/edgelist.txt"),
            0.99,
        )
        .expect("Error while building pedigree: ");

        let model = ABneutral::run(&pedigree, p0uu, p0uu, 1.0, 1000, None).expect("Model failed");
        let result = BootModel::run(&pedigree, &model, p0uu, p0uu, 1.0, 200, None)
            .expect("Bootstrap failed");
        println!("{result}");
        assert_close_10_percent!(model.alpha, 5.7985750419976e-05);
        assert_close_10_percent!(model.beta, 0.00655710970515347);
        assert_close_10_percent!(p0uu, 0.991008120326199);
        assert_close!(model.alpha, 5.7985750419976e-05);
        assert_close!(model.beta, 0.00655710970515347);
        assert_close!(p0uu, 0.991008120326199);

        // assert_close_10_percent!(result.alpha, 1.19025049000535e-05);
        // assert_close_10_percent!(result.beta, 0.00138056339729951);
        // assert_close_10_percent!(result.alpha_beta, 0.594234575518036);
    }
}
