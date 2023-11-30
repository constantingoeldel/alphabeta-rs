use plotters::prelude::*;
use std::path::Path;

pub fn bootstrap(
    alphas: Vec<f64>,
    betas: Vec<f64>,
    output_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let output_file = output_dir.join("bootstrap.png");

    let max = alphas
        .iter()
        .chain(betas.iter())
        .cloned()
        .fold(0.0f64, |a, b| a.max(b));

    let root = BitMapBackend::new(&output_file, (640 * 2, 480 * 2)).into_drawing_area();

    root.fill(&WHITE)?;
    let x_axis = ["Alpha", "Beta"];
    let mut chart = ChartBuilder::on(&root)
        .y_label_area_size(90)
        .x_label_area_size(30)
        .margin(20)
        .caption("Bootstrap Boxplot", ("sans-serif", 40))
        .build_cartesian_2d(x_axis.into_segmented(), 0.0..(max * 1.3) as f32)?;

    chart
        .configure_mesh()
        .y_desc("Epimutation rate")
        .y_labels(10)
        .axis_desc_style(("sans-serif", 30))
        .light_line_style(WHITE)
        .draw()?;

    chart.configure_mesh().draw()?;

    let alpha_quartiles = Quartiles::new(&alphas);
    let beta_quartiles = Quartiles::new(&betas);

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

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

#[cfg(test)]
mod test {
    // use std::path::PathBuf;

    // use super::*;
    // use crate::arguments::Windows;

    // #[test]
    // fn test_metaplot() {
    //     let models: Vec<Analysis> = (0..300).map(|_| Analysis::random()).collect();
    //     let args = Windows::default();

    //     metaplot(&models, &args).unwrap();
    // }

    // #[test]
    // fn test_boxplot() {
    //     let models: Vec<Analysis> = (0..300).map(|_| Analysis::random()).collect();
    //     let alphas = models.iter().map(|m| m.alpha).collect_vec();
    //     let betas = models.iter().map(|m| m.beta).collect_vec();
    //     let out = PathBuf::from("./");
    //     bootstrap(alphas, betas, &out).unwrap();
    // }
}
