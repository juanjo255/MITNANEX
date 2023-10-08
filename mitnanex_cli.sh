#!/bin/bash

## default values
proportion=0.4
threads=4
min_len=-1
max_len=-1
coverage=-1
timestamp=$(date -u +"%Y-%m-%d %T")

# This is for avoiding creating directories
suffix=$(date  "+%Y-%m-%d_%H:%M:%S")
prefix=-1

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
        -r        Prefix name add to every produced file. [input file name]
        -c        Coverage. Minimum coverage per cluster accepted. [-1]
        *         Help.
    "
    exit 1
}

while getopts 'i:t:p:m:M:w:c:r:' opt; do
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
        c)
        coverage=$OPTARG
        ;;
        r)
        prefix=$OPTARG
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

## PREFIX name to use for the resulting files
if [ $prefix -ne -1 ]; then 
    prefix=$(basename $input_file)
    prefix=${prefix%%.*}
fi

## WORKING DIRECTORY
if [ ${wd: -1} = / ]; then 
    wd=$wd"mitnanex_results_"$suffix"/"
else
    wd=$wd"/mitnanex_results_"$suffix"/"
fi

##### FUNCTIONS #####
create_wd(){
## CREATE WORKING DIRECTORY
mkdir $wd
}

subsample(){
### SEQKIT
echo $timestamp': Running seqkit'
seqkit seq -g --threads $threads --min-len $min_len --max-len $max_len \
    $input_file  | \
    seqkit sample --proportion $proportion --threads $threads | \
    seqkit sort --threads $threads --by-length --reverse \
    -o $wd$prefix"_sample.sorted.fastq"
}

mapping(){
### MINIMAP2
echo $timestamp': Running minimap2'
minimap2 -x ava-ont -t $threads --dual=yes --split-prefix $prefix \
    $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq" | \
    fpa keep --containment > $wd$prefix"_containments.paf"
    # fpa drop --dovetail > $wd$prefix"_containments.paf"
}

mt_reads_filt(){
## MITNANEX main
python3 'src/mitnanex.py' $wd$prefix"_sample.sorted.fastq" $wd$prefix"_containments.paf" $coverage
}

first_assembly(){
## FIRST DARFT ASSEMBLY 
flye --scaffold -t $threads --no-alt-contigs --nano-raw $wd$prefix"_putative_mt_reads.fasta" -o $wd"_flye/"
}

contig_selection(){
## SELECT CONTIG AND SUMMON MORE READS
python3 'src/select_contig.py' $wd"_flye/"
}


## START TIMER
start=$SECONDS

### VISAJE INICIAL ###
echo "
    
  __  __ ___ _____ _   _    _    _   _ _______  __
 |  \/  |_ _|_   _| \ | |  / \  | \ | | ____\ \/ /
 | |\/| || |  | | |  \| | / _ \ |  \| |  _|  \  / 
 | |  | || |  | | | |\  |/ ___ \| |\  | |___ /  \ 
 |_|  |_|___| |_| |_| \_/_/   \_\_| \_|_____/_/\_\
                                                  
https://github.com/juanjo255/MITNANEX_PROJECT.git                                     

$timestamp -> Prefix to name resulting files: $prefix

$timestamp -> Working directory: $wd
"

#### PIPELINE ####

create_wd && subsample && mapping #&& mt_reads_filt && first_assembly && contig_selection 

## END TIMER
duration=$(( SECONDS - start ))
echo "Elapsed time: $duration secs."