use alphabeta::*;

fn main() {
    let params = Model::new(0.07);
    let pedigree = Pedigree::from_file("./pedigree.txt");
    println!("{params}");
    let model = ABneutral::run(&pedigree, 0.75, 0.5, 0.7, 1).expect("Model failed");
    let result = BootModel::run(&pedigree, &model, 0.75, 0.5, 0.7, 1).expect("Bootstrap failed");

    println!("##########");
    println!("Model");
    println!("{model}");
    println!("##########");
    println!("Bootstrap");
    println!("{result}");
}
