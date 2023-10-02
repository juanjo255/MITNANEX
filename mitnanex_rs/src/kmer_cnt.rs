use std::cmp;

use itertools::Itertools;

pub fn reverse_compl(seq:&str) -> String {
    let mut rev_comp = String::new();
    for letter in seq.chars(){
        match letter{
            'A' => rev_comp= String::from("T") + &rev_comp,
            'T' => rev_comp= String::from("A") + &rev_comp,
            'C' => rev_comp= String::from("G") + &rev_comp,
            'G' => rev_comp= String::from("C") + &rev_comp,
            _ => panic!("There is some non-canonical nucleotides!")
        };
    };
    return rev_comp
}

pub fn possible_kmers(k: usize, ncl: &str) -> Vec<String> {
    // Create all kmers possible for a value K and keep only canonical kmers

    let total_kmers = (1..=k).map(|_| ncl.chars()).multi_cartesian_product().map(|x| x.iter().join("")).collect_vec();
    let canon_kmers = total_kmers.into_iter().filter(
        |x| {
            &x == cmp::min(&x, &&reverse_compl(&x))
            }
        ).collect_vec();
    return canon_kmers;
}

pub fn count_kmer(k: u8, seq: &str) {
    // Count kmers in a dna sequence using hash table
    todo!()
}
