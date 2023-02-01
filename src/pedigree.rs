use ndarray::{Array2, ArrayView};
pub fn read_pedigree_from_file(filename: &str) -> ndarray::Array2<f64> {
    let mut pedigree = Array2::<f64>::zeros((0, 4));
    let file = std::fs::read_to_string(filename).unwrap();
    file.split('\n').skip(1).for_each(|line| {
        let mut entries = line.split(' ');
        let row = [
            entries.next().unwrap().parse::<f64>().unwrap(),
            entries.next().unwrap().parse::<f64>().unwrap(),
            entries.next().unwrap().parse::<f64>().unwrap(),
            entries.next().unwrap().parse::<f64>().unwrap(),
        ];
        pedigree.push_row(ArrayView::from(&row)).unwrap();
    });
    pedigree
}
