use crate::kmer_cnt::possible_kmers;

pub mod kmer_cnt;

fn main() {
    println!("Mitnanex");
    possible_kmers(3, "ATCG")
}
