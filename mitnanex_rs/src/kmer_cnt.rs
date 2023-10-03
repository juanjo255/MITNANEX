use itertools::Itertools;
use std::{cmp, usize, collections::HashMap, str};
use bio::alphabets::dna;

pub fn reverse_compl(seq: &str) -> String {
    let mut rev_comp = String::new();
    for letter in seq.chars() {
        match letter {
            'A' => rev_comp = String::from("T") + &rev_comp,
            'T' => rev_comp = String::from("A") + &rev_comp,
            'C' => rev_comp = String::from("G") + &rev_comp,
            'G' => rev_comp = String::from("C") + &rev_comp,
            _ => panic!("There is some non-canonical nucleotides!"),
        };
    }
    return rev_comp;
}

pub fn possible_kmers(k: usize, ncl: &str) -> Vec<String> {
    // Create all kmers possible for a value K and keep only canonical kmers

    let total_kmers:Vec<String> = (1..=k)
        .map(|_| ncl.chars())
        .multi_cartesian_product()
        .map(|x| x.iter().join(""))
        .collect_vec();
    let canon_kmers:Vec<String> = total_kmers
        .into_iter()
        .filter(|x| &x == cmp::min(&x, &&reverse_compl(&x)))
        .collect_vec();
    return canon_kmers;
}

pub fn count_kmer(k: usize, seq: &str) -> HashMap<String, f32> {
    // Count kmers in a dna sequence using hash table
    let len_seq:usize = seq.len();
    let num_kmer_expected:usize = len_seq - k + 1;
    let possible_kmers_list:Vec<String> = possible_kmers(k, "ATCG");
    let mut kmer_counter:HashMap<String, f32> = HashMap::from_iter(possible_kmers_list.into_iter().map(|x| (x, 1.0)));

    for kmer in 0..num_kmer_expected {
        let kmer_for = &seq[kmer..kmer+k];
        let kmer_rev = dna::revcomp(kmer_for.as_bytes());
        let kmer_rev = str::from_utf8(&kmer_rev).unwrap();
        let canon_kmer = cmp::min(kmer_for, &kmer_rev);
        if let Some(x) = kmer_counter.get_mut(canon_kmer){
            *x+= 1.0;
        }else {
            panic!("There is non-canonical nucleotides!")
        };
    };

    // Normalize kmer count by read length
    for kmer in kmer_counter.clone().keys(){
        if let Some(x) = kmer_counter.get_mut(kmer){
            *x = *x / len_seq as f32;
        };
    };
    return kmer_counter
}

pub fn get_kmer_profiles (reads:&Vec<&str>, reads_file:&str) -> Vec<f32>{
    
    todo!()
}