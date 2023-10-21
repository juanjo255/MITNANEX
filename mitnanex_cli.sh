#!/bin/bash
## default values
proportion=0.4
threads=4
min_len=-1
max_len=-1
coverage=-1
timestamp=$(date -u +"%Y-%m-%d %T")
map_identity=0.6
min_qual=-1
wd="./"
output_dir='mitnanex_results/'

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
        -m        Min-len. Filter reads by minimun length. Read seqkit seq documentation. [-1].
        -M        Max-len. Filter reads by maximun length. Read seqkit seq documentation. [-1].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./mitnanex_results].
        -r        Prefix name add to every produced file. [input file name].
        -c        Coverage. Minimum coverage per cluster accepted. [-1].
        -d        Different output directory. Create a different output directory every run (it uses the date and time).
        -s        Mapping identity. Minimun identity between two reads to be store in the same cluster.[0.6]
        -q        Min mapping quality (>=). This is for samtools. [-1].
        *         Help.
    "
    exit 1
}

while getopts 'i:t:p:m:M:w:c:r:s:q:d' opt; do
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
        s)
        map_identity=$OPTARG
        ;;
        q)
        min_qual=$OPTARG
        ;;
        d)
        output_dir="mitnanex_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        *)
        mitnanex_help
        ;;
    esac 
done

# Check if required arguments are provided
if [ -z "$input_file" ];
then
  echo "Error: Input file is required."
  mitnanex_help
fi

## PREFIX name to use for the resulting files
if [ -z $prefix ];
then 
    prefix=$(basename $input_file)
    prefix=${prefix%%.*}
fi

if [ ${wd: -1} = / ];
then 
    wd=$wd$output_dir
else
    wd=$wd"/"$output_dir
fi

##### FUNCTIONS #####
create_wd(){
## CREATE WORKING DIRECTORY
if [ -d $wd ]
then
  echo "Rewriting directory..."
  echo " "
else 
    mkdir $wd
fi
}

subsample(){
### SEQKIT
echo " "
echo $timestamp': Step 1: Sampling with seqkit'
echo " "
seqkit seq -g --threads $threads --min-len $min_len --max-len $max_len $input_file | \
    seqkit sample --proportion $proportion --threads $threads -o $wd$prefix"_sample.sorted.fastq"
    
}

trim_adapters(){
## $1 input
## $2 output
    echo " "
    echo $timestamp': Step 2: Trimming adapters with porechop'
    echo " "
    porechop --verbosity 0 -t $threads -i $1 -o $2 
}

sort_file(){
    echo " "
    echo $timestamp': Step 3: Sorting file with seqkit'
    echo " "
    seqkit sort --threads $threads --by-length --reverse -o $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq"
}

reads_overlap(){
### MINIMAP2
echo " "
echo $timestamp': Step 4:  Running minimap2'
echo " "
minimap2 -x ava-ont -t $threads --dual=yes --split-prefix $prefix \
    $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq" | \
    fpa drop --internalmatch --length-lower $min_len > $wd$prefix"_containments.paf"
}

mt_reads_filt(){
## MITNANEX main
echo " "
echo $timestamp': Step 5: Running MITNANEX'
echo " "
python3 main.py $wd$prefix"_sample.sorted.fastq" $wd$prefix"_containments.paf" $coverage $map_identity $wd$prefix"_putative_mt_reads.fasta"
}

first_assembly(){
## FIRST DARFT ASSEMBLY 
echo " "
echo $timestamp': Step 7: Running Miniasm'
echo " "
minimap2 -x ava-ont -t $threads --dual=yes --split-prefix $prefix \
    $wd$prefix"_putative_mt_reads.fasta" $wd$prefix"_putative_mt_reads.fasta" | \
    miniasm -f $wd$prefix"_putative_mt_reads.fasta" - > $wd$prefix"_first_draft_asm.gfa"
}

statistics(){
## STATISTICS OF FIRST DRAFT ASSEMBLY
    echo " "  
    echo $timestamp': Step 8: Computing statistics with gfastats'
    echo " "
    echo "#### DESCRIPTION OF CONTIGS FOUND #### "
    gfastats  --seq-report --discover-paths $wd$prefix"_first_draft_asm.gfa"
    echo " "
}

gfa2fasta(){
## SELECT CONTIG AND SUMMON MORE READS
    ## python3 src/select_contig.py $wd$prefix"_flye/" $wd$prefix"_first_draft_mt_assembly.fasta"
    echo " "
    echo $timestamp': Step 9: Converting gfa to fasta'
    echo " "
    ### convert form gfa to fasta
    gfastats --discover-paths $wd$prefix"_first_draft_asm.gfa" -o $wd$prefix"_first_draft_asm.fasta" > /dev/null
}

collecting_mt_reads(){
## $1 reference minimap
## $2 input file minimap
## $3 output minimap
## $4 output samtools

## USING MINIASM ASSEMBLY COLLECT MORE READS
    echo " "
    echo $timestamp': Step 10: Collecting more reads using the draft assembly with minimap2'
    echo " "
    ### Map reads to the unitig formed by miniasm
    minimap2 -ax map-ont --split-prefix $prefix  $1 $2 -o $3
    ### Get correctly mapped reads
    echo " "
    echo $timestamp" : Reads collected for final assembly"
    samtools view --min-MQ $min_qual -F4 --bam $3 | samtools fastq > $4
    echo " "
}

final_assembly(){
    echo $timestamp': Step 11: Running final assembly with Flye'
    flye --scaffold -t $threads --iterations 5 --no-alt-contigs --nano-raw $wd$prefix"_collected_reads.fastq" -o $wd$prefix"_flye/"
}



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

## START TIMER
start=$SECONDS

#### PIPELINE ####
create_wd && subsample && trim_adapters $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq" \
&& sort_file && reads_overlap && mt_reads_filt && first_assembly && statistics && gfa2fasta \
#&& collecting_mt_reads $wd$prefix"_first_draft_asm.fasta" $input_file $wd$prefix"_align.sam" $wd$prefix"_collected_reads.fastq" \
#&& final_assembly

#mt_reads_filt && first_assembly && statistics && gfa2fasta \ 
#&& collecting_mt_reads $wd$prefix"_first_draft_asm.fasta" $input_file $wd$prefix"_align.sam" $wd$prefix"_collected_reads.fastq" \
#&& final_assembly

echo ""
echo "### MITNANEX finished ###"
echo ""
echo "Final assembly is in" $wd$prefix"_flye""/assembly.fasta"
echo "Putative mitochondrial reads are in "$wd$prefix"_collected_reads.fastq"

## END TIMER
duration=$(( SECONDS - start ))
echo "" 
echo "$timestamp -> Elapsed time: $duration secs."
echo "" 