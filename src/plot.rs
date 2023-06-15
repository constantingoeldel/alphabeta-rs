use std::path::PathBuf;

use crate::{
    arguments::Windows,
    structs::{Model, ModelWithSD},
    *,
};
use itertools::Itertools;
use plotters::prelude::*;

pub fn metaplot(windows: &[ModelWithSD], args: &Windows) -> Result<()> {
    let output_file = args.output_dir.join("metaplot.png");
    dbg!(&output_file);

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
            windows
                .iter()
                .enumerate()
                .map(|(i, w)| (i as f32, w.alpha as f32)),
            &RED,
        ))?
        .label("Alpha")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));

    chart
        .draw_series(LineSeries::new(
            windows
                .iter()
                .enumerate()
                .map(|(i, w)| (i as f32, w.beta as f32)),
            &BLUE,
        ))?
        .label("Beta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

    let confidence_interval_alpha: Vec<(f32, f32)> = windows
        .iter()
        .enumerate()
        .map(|(i, w)| {
            vec![
                (i as f32, (w.alpha - w.sd_alpha) as f32),
                (i as f32, (w.alpha + w.sd_alpha) as f32),
            ]
        })
        .concat();
    chart.draw_series(std::iter::once(Polygon::new(
        confidence_interval_alpha,
        RED.mix(0.2),
    )))?;

    let confidence_interval_beta: Vec<(f32, f32)> = windows
        .iter()
        .enumerate()
        .map(|(i, w)| {
            vec![
                (i as f32, (w.beta - w.sd_beta) as f32),
                (i as f32, (w.beta + w.sd_beta) as f32),
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

pub fn bootstrap(alphas: &[f64], betas: &[f64]) -> Result<()> {
    let output_file = PathBuf::from("bootstrap.png");
    dbg!(&output_file);

    let max = alphas
        .iter()
        .chain(betas.iter())
        .cloned()
        .fold(0.0f64, |a, b| a.max(b));

    let root = BitMapBackend::new(&output_file, (640 * 2, 480 * 2)).into_drawing_area();

    root.fill(&WHITE)?;
    let x_axis = ["Alpha", "Beta"];
    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(40)
        .y_label_area_size(80)
        .caption("Bootstrap Boxplot", ("sans-serif", 20))
        .build_cartesian_2d(x_axis.into_segmented(), 0.0..(max * 1.3) as f32)?;

    chart
        .configure_mesh()
        .x_desc("Ping, ms")
        .y_desc("Host")
        .y_labels(10)
        .light_line_style(WHITE)
        .draw()?;

    chart.configure_mesh().draw()?;

    let alpha_quartiles = Quartiles::new(alphas);
    let beta_quartiles = Quartiles::new(betas);

    chart.draw_series(vec![
        Boxplot::new_vertical(SegmentValue::CenterOf(&"Alpha"), &alpha_quartiles)
            .width(100)
            .style(RED)
            .whisker_width(1.0),
        Boxplot::new_vertical(SegmentValue::CenterOf(&"Beta"), &beta_quartiles)
            .width(100)
            .style(BLUE)
            .whisker_width(1.0),
    ])?;

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::arguments::Windows;

    #[test]
    fn test_metaplot() {
        let models: Vec<ModelWithSD> = (0..300).map(|_| ModelWithSD::random()).collect();
        let args = Windows::default();

        metaplot(&models, &args).unwrap();
    }

    #[test]
    fn test_boxplot() {
        let models: Vec<ModelWithSD> = (0..300).map(|_| ModelWithSD::random()).collect();
        let alphas = models.iter().map(|m| m.alpha).collect_vec();
        let betas = models.iter().map(|m| m.beta).collect_vec();

        bootstrap(&alphas, &betas).unwrap();
    }
}
