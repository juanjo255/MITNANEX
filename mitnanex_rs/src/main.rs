use crate::kmer_cnt::count_kmer;

pub mod kmer_cnt;

fn main() {
    println!("Mitnanex");
    count_kmer(3, "ACTTCA");
}
