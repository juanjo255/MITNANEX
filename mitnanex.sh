#!/bin/bash

## default values
proportion=0.4
threads=4

## Help message
mitnanex_help() {
    echo "
        MITNANEX
        Version: 1.0
        https://github.com/juanjo255/MITNANEX_PROJECT.git

        Usage: $mitnanex [options] FASTQ

        Flags:
        -t        Threads
        -p        Proportion. for sampling with seqkit. for information read seqkit sample documentation
        *         Help
    "
    exit 1
}

while getopts 'i:t:p' opt; do
    case $opt in
        i)
        input_file=$OPTARG
        ;;
        t)
        threads=$OPTARG
        ;;
        p)
        proportion=$OPTARG
        ;;
        *)
        mitnanex_help
        ;;
    esac 


## prefix name to use for the resulting files 
prefix=${input_file%%.*}

seqkit sample --proportion $proportion --threads $threads \
    $input_file  \ 
    -o $prefix"_sample.fastq" 

# seqkit sort --threads $threads --by-length --reverse \ 
#     s_cervisae_CEN.PK113-7D_SRR5892449_reads_sample.fastq \ 
#     -o $prefix"_sample.sorted.fastq"

# minimap2 -x ava-ont -t $threads --dual=yes --split-prefix \ 
#     temp_name $prefix"_sample.sorted.fastq" $prefix"_sample.sorted.fastq" | \ 
#     fpa keep --containment > $prefix"_containments.paf"
