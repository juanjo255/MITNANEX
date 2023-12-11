use bio::alphabets::dna;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use indexmap::IndexMap;
use itertools::Itertools;
use pyo3::prelude::*;
use std::fs::File;
use std::{cmp, str, usize};

pub fn possible_kmers(k: usize, ncl: &str) -> Vec<String> {
    // Create all kmers possible for a value K and keep only canonical kmers

    let total_kmers: Vec<String> = (1..=k)
        .map(|_| ncl.chars())
        .multi_cartesian_product()
        .map(|x| x.iter().join(""))
        .collect_vec();
    let canon_kmers: Vec<String> = total_kmers
        .into_iter()
        .filter(|x| {
            &x == cmp::min(
                &x,
                &&str::from_utf8(&dna::revcomp(x.as_bytes()))
                    .unwrap()
                    .to_string(),
            )
        })
        .collect_vec();
    return canon_kmers;
}

pub fn count_kmer(k: usize, seq: &str) -> IndexMap<String, f32> {
    // Count kmers in a dna sequence using hash table
    let len_seq: usize = seq.len();
    let num_kmer_expected: usize = len_seq - k + 1;
    let possible_kmers_list: Vec<String> = possible_kmers(k, "ATCG");
    let mut kmer_counter: IndexMap<String, f32> =
        IndexMap::from_iter(possible_kmers_list.into_iter().map(|x| (x, 1.0)));

    for kmer in 0..num_kmer_expected {
        let kmer_for = &seq[kmer..kmer + k];
        let kmer_rev = dna::revcomp(kmer_for.as_bytes());
        let kmer_rev = str::from_utf8(&kmer_rev).unwrap();
        let canon_kmer = cmp::min(kmer_for, &kmer_rev);
        if let Some(x) = kmer_counter.get_mut(canon_kmer) {
            *x += 1.0;
        } else {
            panic!("There is non-canonical nucleotides!")
        };
    }

    // Normalize kmer count by read length
    for kmer in kmer_counter.clone().keys() {
        if let Some(x) = kmer_counter.get_mut(kmer) {
            *x = *x / len_seq as f32;
        };
    }
    return kmer_counter;
}
//reads_ids:&Vec<&str>,

#[pyfunction]
pub fn get_kmer_profiles(reads_ids: Vec<&str>, reads_file: &str, k: usize) -> (Vec<Vec<f32>>, Vec<String>) {
    let mut reader = fastq::Reader::new(
        File::open(reads_file)
            .expect("Ups! I could not open the file. Check the path and the file"),
    );
    let mut record = fastq::Record::new();
    reader.read(&mut record).expect("Failed to parse record");
    let mut kmer_profiles: Vec<Vec<f32>> = Vec::new();
    let mut kmer_profiles_ids: Vec<String> = Vec::new();

    // Go through each record in the fastQ
    while !record.is_empty() {
        let check = record.check();
        if check.is_err() {
            panic!("I got a bad record!")
        }
        if reads_ids.contains(&record.id()) {
            // obtain sequence
            let seq = record.seq();
            let seq = str::from_utf8(&seq).unwrap();
            let kmer_comp = count_kmer(k, seq);
            kmer_profiles.push(kmer_comp.into_values().collect_vec());
            kmer_profiles_ids.push(record.id().to_string());
        }
        // Update record
        // Equivalent to readlines in python
        reader.read(&mut record).expect("Failed to parse record");
    }
    return (kmer_profiles, kmer_profiles_ids);
}

/// A Python module implemented in Rust.
#[pymodule]
fn utils_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_kmer_profiles, m)?)?;
    Ok(())
}
