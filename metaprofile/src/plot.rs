use alphabeta::Analysis;
use itertools::Itertools;
use plotters::prelude::*;

use crate::{Metaprofile, Return};

pub fn metaplot(
    analyses: &[Analysis],
    args: &Metaprofile,
) -> Result<(), Box<dyn std::error::Error>> {
    let output_file = args.output_dir.join("metaplot.png");

    let root = BitMapBackend::new(&output_file, (640 * 2, 480 * 2)).into_drawing_area();

    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Metaplot", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f32..300f32, 0f32..0.01f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            analyses
                .iter()
                .enumerate()
                .map(|(i, w)| (i as f32, w.alpha as f32)),
            &RED,
        ))?
        .label("Alpha")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));

    chart
        .draw_series(LineSeries::new(
            analyses
                .iter()
                .enumerate()
                .map(|(i, w)| (i as f32, w.beta as f32)),
            &BLUE,
        ))?
        .label("Beta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

    let confidence_interval_alpha: Vec<(f32, f32)> = analyses
        .iter()
        .enumerate()
        .map(|(i, w)| {
            vec![
                (i as f32, (w.ci_alpha.1) as f32),
                (i as f32, (w.ci_alpha.0) as f32),
            ]
        })
        .concat();
    chart.draw_series(std::iter::once(Polygon::new(
        confidence_interval_alpha,
        RED.mix(0.2),
    )))?;

    let confidence_interval_beta: Vec<(f32, f32)> = analyses
        .iter()
        .enumerate()
        .map(|(i, w)| {
            vec![
                (i as f32, (w.ci_beta.1) as f32),
                (i as f32, (w.ci_beta.0) as f32),
            ]
        })
        .concat();
    chart.draw_series(std::iter::once(Polygon::new(
        confidence_interval_beta,
        BLUE.mix(0.2),
    )))?;

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}
