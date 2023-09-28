#!/bin/bash

## default values
proportion=0.4
threads=4
min_len=500
timestamp=$(date -u +"%Y-%m-%d %T")

## Help message
mitnanex_help() {
    echo "
    MITNANEX
    Version: 1.0
    https://github.com/juanjo255/MITNANEX_PROJECT.git

    Usage: mitnanex.sh [options] FASTQ

    Options:
        -i        Input file.
        -t        Threads.
        -p        Proportion. For sampling with seqkit. Read seqkit sample documentation.
        -m        Min-len. Filter reads by minimun length. Read seqkit seq documentation.
        *         Help.
    "
    exit 1
}

while getopts 'i:t:p:m' opt; do
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
        min_len=$OPTARG
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

echo '
    
  __  __ ___ _____ _   _    _    _   _ _______  __
 |  \/  |_ _|_   _| \ | |  / \  | \ | | ____\ \/ /
 | |\/| || |  | | |  \| | / _ \ |  \| |  _|  \  / 
 | |  | || |  | | | |\  |/ ___ \| |\  | |___ /  \ 
 |_|  |_|___| |_| |_| \_/_/   \_\_| \_|_____/_/\_\
                                                  
https://github.com/juanjo255/MITNANEX_PROJECT.git                                     
'
                                                              

### SEQKIT
echo $timestamp': Running seqkit'
seqkit seq -g --threads $threads --min-len $min_len \
    $input_file  | \
    seqkit sample --proportion $proportion --threads $threads | \
    seqkit sort --threads $threads --by-length --reverse \
    -o $prefix"_sample.sorted.fastq"

### MINIMAP2
echo $timestamp': Running minimap2'
minimap2 -x ava-ont -t $threads --dual=yes --split-prefix $prefix \
    $prefix"_sample.sorted.fastq" $prefix"_sample.sorted.fastq" | \
    fpa keep --containment > $prefix"_containments.paf"

