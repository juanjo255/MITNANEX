use itertools::Itertools;

pub fn possible_kmers (k:u8, ncl:&str){
    let total_kmers = (1..=k).cartesian_product(ncl.chars());
    println!("{:?}",total_kmers);
}

pub fn count_kmer (k:u8,seq:&str) {
    // Count kmers in a dna sequence using hash table
    todo!()
}