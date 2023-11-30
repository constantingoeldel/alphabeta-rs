use std::{path::Path, fs::File, fmt::{Display, Debug}, ops::Div, io::Write};

use ndarray::{Array2, ArrayView1, array};
use ndarray_stats::{interpolate::Linear, Quantile1dExt};
use noisy_float::types::n64;

/// 2D-Array containing the results of all the iterations of the bootstrapping analysis.
/// 
/// Dimensions: 7 x n_bootze
/// 
/// Columns: Alpha, Beta, Weight, Intercept, Pr(MM), Pr(UM), Pr(UU)
pub struct RawAnalysis(pub Array2<f64>);


#[derive(Debug, Clone)]
pub struct Analysis {
    pub alpha: f64,
    pub beta: f64,
    pub alphabeta: f64,
    pub weight: f64,
    pub intercept: f64,

    pub pr_mm: f64,
    pub pr_um: f64,
    pub pr_uu: f64,

    pub sd_alpha: f64,
    pub sd_beta: f64,
    pub sd_alphabeta: f64,
    pub sd_weight: f64,
    pub sd_intercept: f64,

    pub sd_pr_mm: f64,
    pub sd_pr_um: f64,
    pub sd_pr_uu: f64,

    pub ci_alpha: CI,
    pub ci_beta: CI,
    /// Beta/Alpha
    pub ci_alphabeta: CI,
    pub ci_weight: CI,
    pub ci_intercept: CI,

    pub ci_pr_mm: CI,
    pub ci_pr_um: CI,
    pub ci_pr_uu: CI,
}

impl RawAnalysis{
  pub  fn analyze(&self) -> Analysis {
        
    
    let results = self.0.view();
    let alphabeta = results.column(1).div(&results.column(0));

    let ci = |val: ArrayView1<f64>| {
        let q = val
            .map(|x| n64(*x))
            .to_owned()
            .quantiles_mut(&array![n64(0.025), n64(0.975)], &Linear)
            .unwrap();

        CI(q[0].raw(), q[1].raw())
    };

   Analysis {
        alpha: results.column(0).mean().unwrap(),
        beta: results.column(1).mean().unwrap(),
        alphabeta: alphabeta.mean().unwrap(),
        weight: results.column(2).mean().unwrap(),
        intercept: results.column(3).mean().unwrap(),
        pr_mm: results.column(4).mean().unwrap(),
        pr_um: results.column(5).mean().unwrap(),
        pr_uu: results.column(6).mean().unwrap(),

        sd_alpha: results.column(0).std(1.0),
        sd_beta: results.column(1).std(1.0),

        sd_alphabeta: alphabeta.std(1.0),

        sd_weight: results.column(2).std(1.0),
        sd_intercept: results.column(3).std(1.0),

        sd_pr_mm: results.column(4).std(1.0),
        sd_pr_um: results.column(5).std(1.0),
        sd_pr_uu: results.column(6).std(1.0),

        ci_alpha: ci(results.column(0)),
        ci_beta: ci(results.column(1)),
        ci_alphabeta: ci(alphabeta.view()),
        ci_weight: ci(results.column(2)),
        ci_intercept: ci(results.column(3)),

        ci_pr_mm: ci(results.column(4)),
        ci_pr_um: ci(results.column(5)),
        ci_pr_uu: ci(results.column(6)),
    }
    }
}

impl Analysis {
    pub fn to_file<P>(&self, path: P) -> std::io::Result<()> where P: AsRef<Path> + Debug {
        println!("Writing model to file: {:?}", path);
        let mut file = File::create(path).unwrap();
        // TODO: Use the display implementation from below
        let  content = format!(
            "Alpha\t{}\nBeta\t{}\nAlphaBeta\t{}\nWeight\t{}\nIntercept\t{}\nPrMM\t{}\nPrUM\t{}\nPrUU\t{}\nSDAlpha\t{}\nSDBeta\t{}\nSDAlphaBeta\t{}\nSDWeight\t{}\nSDIntercept\t{}\nSDPrMM\t{}\nSDPrUM\t{}\nSDPrUU\t{}\nCIAlpha\t{}-{}\nCIBeta\t{}-{}\nCIAlphaBeta\t{}-{}\nCIWeight\t{}-{}\nCIIntercept\t{}-{}\nCIPrMM\t{}-{}\nCIPrUM\t{}-{}\nCIPrUU\t{}-{}\n",
            self.alpha,
            self.beta,
            self.alphabeta,
            self.weight,
            self.intercept,
            self.pr_mm,
            self.pr_um,
            self.pr_uu,
            self.sd_alpha,
            self.sd_beta,
            self.sd_alphabeta,
            self.sd_weight,
            self.sd_intercept,
            self.sd_pr_mm,
            self.sd_pr_um,
            self.sd_pr_uu,
            self.ci_alpha.0,
            self.ci_alpha.1,
            self.ci_beta.0,
            self.ci_beta.1,
            self.ci_alphabeta.0,
            self.ci_alphabeta.1,
            self.ci_weight.0,
            self.ci_weight.1,
            self.ci_intercept.0,
            self.ci_intercept.1,
            self.ci_pr_mm.0,
            self.ci_pr_mm.1,
            self.ci_pr_um.0,
            self.ci_pr_um.1,
            self.ci_pr_uu.0,
            self.ci_pr_uu.1,
            

        );

        file.write_all(content.as_bytes())
    }
}

impl Display for Analysis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, 
            "Alpha\t{}\nBeta\t{}\nAlphaBeta\t{}\nWeight\t{}\nIntercept\t{}\nPrMM\t{}\nPrUM\t{}\nPrUU\t{}\nSDAlpha\t{}\nSDBeta\t{}\nSDAlphaBeta\t{}\nSDWeight\t{}\nSDIntercept\t{}\nSDPrMM\t{}\nSDPrUM\t{}\nSDPrUU\t{}\nCIAlpha\t{}-{}\nCIBeta\t{}-{}\nCIAlphaBeta\t{}-{}\nCIWeight\t{}-{}\nCIIntercept\t{}-{}\nCIPrMM\t{}-{}\nCIPrUM\t{}-{}\nCIPrUU\t{}-{}\n",
            self.alpha,
            self.beta,
            self.alphabeta,
            self.weight,
            self.intercept,
            self.pr_mm,
            self.pr_um,
            self.pr_uu,
            self.sd_alpha,
            self.sd_beta,
            self.sd_alphabeta,
            self.sd_weight,
            self.sd_intercept,
            self.sd_pr_mm,
            self.sd_pr_um,
            self.sd_pr_uu,
            self.ci_alpha.0,
            self.ci_alpha.1,
            self.ci_beta.0,
            self.ci_beta.1,
            self.ci_alphabeta.0,
            self.ci_alphabeta.1,
            self.ci_weight.0,
            self.ci_weight.1,
            self.ci_intercept.0,
            self.ci_intercept.1,
            self.ci_pr_mm.0,
            self.ci_pr_mm.1,
            self.ci_pr_um.0,
            self.ci_pr_um.1,
            self.ci_pr_uu.0,
            self.ci_pr_uu.1,
            

        )
    }
}

/// Lower and upper bounds of the confidence interval for a given probability
///
/// (lower, upper)
///
/// (0.025, 0.975) for 95% CI
#[derive(Debug, Clone)]
pub struct CI(pub f64, pub f64);

