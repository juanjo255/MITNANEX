#!/bin/bash

## START TIMER
start=$SECONDS

## default values
proportion=0.4
threads=4
min_len=-1
max_len=-1
timestamp=$(date -u +"%Y-%m-%d %T")
wd="./mitnanex_results"

## Help message
mitnanex_help() {
    echo "
    MITNANEX
    Version: 1.0
    https://github.com/juanjo255/MITNANEX_PROJECT.git

    Usage: mitnanex.sh [options] FASTQ

    Options:
        -i        Input file.
        -t        Threads. [4].
        -p        Proportion. For sampling with seqkit. Read seqkit sample documentation. [0.4].
        -m        Min-len. Filter reads by minimun length. Read seqkit seq documentation. [500].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./mitnanex_results]
        *         Help.
    "
    exit 1
}

while getopts 'i:t:p:m:M:w:' opt; do
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
        m)
        min_len=$((OPTARG))
        ;;
        M)
        max_len=$((OPTARG))
        ;;
        w)
        wd=$OPTARG
        ;;
        *)
        mitnanex_help
        ;;
    esac 
done

# Check if required arguments are provided
if [ -z "$input_file" ]; then
  echo "Error: Input file is required."
  mitnanex_help
fi

## prefix name to use for the resulting files 
prefix=${input_file%%.*}


## pipeline body

echo "
    
  __  __ ___ _____ _   _    _    _   _ _______  __
 |  \/  |_ _|_   _| \ | |  / \  | \ | | ____\ \/ /
 | |\/| || |  | | |  \| | / _ \ |  \| |  _|  \  / 
 | |  | || |  | | | |\  |/ ___ \| |\  | |___ /  \ 
 |_|  |_|___| |_| |_| \_/_/   \_\_| \_|_____/_/\_\
                                                  
https://github.com/juanjo255/MITNANEX_PROJECT.git                                     

Prefix to name resulting files: $prefix
Working directory: $wd
"


## CREATE WORKING DIRECTORY
mkdir $wd &&

### SEQKIT
echo $timestamp': Running seqkit'
seqkit seq -g --threads $threads --min-len $min_len --max-len $max_len \
    $input_file  | \
    seqkit sample --proportion $proportion --threads $threads | \
    seqkit sort --threads $threads --by-length --reverse \
    -o $wd$prefix"_sample.sorted.fastq" &&

### MINIMAP2
echo $timestamp': Running minimap2'
minimap2 -x ava-ont -t $threads --dual=yes --split-prefix $prefix \
    $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq" | \
    fpa keep --containment > $wd$prefix"_containments.paf" &&
    # fpa drop --dovetail > $wd$prefix"_containments.paf"

python3 mitnanex.py &&

#flye --scaffold -t 5 --nano-raw $wd$prefix"_putative_mt_reads.fasta" -o $wd_flye 


## END TIMER
duration=$(( SECONDS - start ))

echo "Elapsed time: $duration"