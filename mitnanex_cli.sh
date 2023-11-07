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
flye_mode='--nano-hq'
output_dir='mitnanex_results/'

## Help message
mitnanex_help() {
    echo "
    MITNANEX - MITochondrial NANopore reads EXtractor

    Author:
    Juan Picon Cossio

    Version: 0.1

    Usage: mitnanex.sh [options] FASTQ

    Options:
        -i        Input file. [required]
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
        -f        Flye mode. [--nano-hq]
        -g        GenomeSize. This is your best estimation of the mitogenome for read correction with Canu. [required]
        *         Help.
    "
    exit 1
}

while getopts 'i:t:p:m:M:w:c:r:s:q:f:g:d' opt; do
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
        f)
        flye_mode=$OPTARG
        ;;
        g)
        genomeSize=$OPTARG
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

# Check if required arguments are provided
if [ -z "$genomeSize" ];
then
  echo "Error: genome size is required."
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
    echo $timestamp": Rewriting directory..."
    echo " "
    else 
        echo $timestamp": Creating directory..."
        echo " "
        mkdir $wd
    fi
}

subsample(){
## SUBSAMPLE
    echo " "
    echo $timestamp': Step 1: Sampling with seqkit'
    echo " "
    seqkit seq -g --threads $threads --min-len $min_len --max-len $max_len $input_file | \
        seqkit sample --proportion $proportion --threads $threads -o $wd$prefix"_sample.sorted.fastq"
    
}

trim_adapters(){
## $1 input
## $2 output
## TRIM ADAPTERS
    echo " "
    echo $timestamp': Trimming adapters with porechop'
    echo " "
    porechop --verbosity 0 -t $threads -i $1 -o $2 
}

sort_file(){
## SORT FILE
    echo " "
    echo $timestamp': Sorting file with seqkit'
    echo " "
    seqkit sort --threads $threads --by-length --reverse \
    -o $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq"
}

reads_overlap(){
## FIND OVERLAP BETWEEN READS
## --split-prefix $prefix
    echo " "
    echo $timestamp': Looking for overlaps with minimap2'
    echo " "
    minimap2 -x ava-ont -t $threads --dual=yes \
    $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq" | \
        fpa drop --internalmatch --length-lower $min_len > $wd$prefix".paf"
}

mt_reads_filt(){
## MITNANEX
    echo " "
    echo $timestamp': Clustering and discriminating potential mt reads with MITNANEX'
    echo " "
    python3 main.py $wd$prefix"_sample.sorted.fastq" $wd$prefix".paf" $coverage $map_identity $wd$prefix"_putative_mt_reads.fasta"
}

first_assembly(){
## ASSEMBLE WITH MINIASM
    echo " "
    echo $timestamp': Generating draft unitigs with Miniasm'
    echo " "
    minimap2 -x ava-ont -t $threads --dual=yes \
    $wd$prefix"_putative_mt_reads.fasta" $wd$prefix"_putative_mt_reads.fasta" | \
        miniasm -S6 -f $wd$prefix"_putative_mt_reads.fasta" - > $wd$prefix"_first_draft_asm.gfa"
    #> $wd$prefix"_miniasm.paf"
    #miniasm -S6 -f $wd$prefix"_putative_mt_reads.fasta" $wd$prefix"_miniasm.paf" > $wd$prefix"_first_draft_asm.gfa"
}

gfa2fasta(){
## Convert from gfa and fasta and compute descriptive statistics
    ## STATISTICS OF FIRST DRAFT ASSEMBLY
    echo " "  
    echo $timestamp': Computing statistics with gfastats and generating fasta'
    echo " "
    echo "#### DESCRIPTION OF CONTIGS FOUND #### "
    gfastats  --seq-report --discover-paths $wd$prefix"_first_draft_asm.gfa"
    echo " "
    ### convert form gfa to fasta
    gfastats --discover-paths $wd$prefix"_first_draft_asm.gfa" -o $wd$prefix"_first_draft_asm.fasta" > /dev/null
}

collecting_mt_reads(){
## USING MINIASM ASSEMBLY COLLECT MORE READS
    ## $1 reference minimap
    ## $2 input file minimap
    ## $3 output minimap
    ## $4 output samtools
    echo " "
    echo $timestamp': Recruiting more reads with minimap2'
    echo " "
    ### Map reads to the unitig formed by miniasm
    minimap2 -ax map-ont --split-prefix $prefix --secondary=no  $1 $2 | samtools view -b --min-MQ $min_qual -F4 -T $1 -o $3
    ### Get correctly mapped reads
    echo " "
    echo $timestamp" : Reads collected"
    samtools fastq $3 > $4
    echo " "
}

second_assembly(){
## ASSEMBLE WITH FLYE
    echo $timestamp': Generating a draft assembly with Flye'
    flye -t $threads --meta \
        $flye_mode $wd$prefix"_collected_reads.filtlong.fastq" -o $wd$prefix"_flye_asm/"
}

select_contig(){
## SELECT A CONTIG GENERATED BY FLYE
    echo $timestamp': Selecting contig with greatest coverage and length'
    python3 src/select_contig.py $wd$prefix"_flye_asm/" $wd$prefix"_second_draft_asm.fasta"
}

correct_reads(){
## CORRECT COLLECTED READS WITH CANU
    echo ""
    echo $timestamp': Correcting reads with Canu'
    mkdir $wd"canu_correction/"
    canu -correct -d $wd"canu_correction/" -p $prefix"_collected_reads" genomeSize=$genomeSize \
        -nanopore $wd$prefix"_collected_reads.filtlong.fastq" 2> /dev/null \
        && mv $wd"canu_correction/"$prefix"_collected_reads.correctedReads.fasta.gz" $wd$prefix"_collected_reads.correctedReads.fasta.gz" \
        && gunzip $wd$prefix"_collected_reads.correctedReads.fasta"
}

polish_asm(){
## POLISHED THE FLYE ASM USING FLYE POLISHING TOOL

    echo ""
    echo $timestamp': Polishing assembly with Canu corrected reads'
    flye -t $threads $flye_mode $wd$prefix"_collected_reads.correctedReads.fasta"  \
        --polish-target $wd$prefix"_second_draft_asm.fasta"  \
        --iterations 3 -o $wd$prefix"flye_asm_polish_canuCorrected"

    mv $wd$prefix"flye_asm_polish_canuCorrected/polished_3.fasta" $wd$prefix"_final_draft_asm.fasta"
}

quality_control(){
## QUALITY CONTROL USING FILTLONG
    echo ""
    echo $timestamp': Tossing up bad reads with filtlong'
    filtlong --min_length $min_len --keep_percent 90 $wd$prefix"_collected_reads.fastq" > $wd$prefix"_collected_reads.filtlong.fastq"
}

### VISAJE INICIAL ###
echo "
    
  __  __ ___ _____ _   _    _    _   _ _______  __
 |  \/  |_ _|_   _| \ | |  / \  | \ | | ____\ \/ /
 | |\/| || |  | | |  \| | / _ \ |  \| |  _|  \  / 
 | |  | || |  | | | |\  |/ ___ \| |\  | |___ /  \ 
 |_|  |_|___| |_| |_| \_/_/   \_\_| \_|_____/_/\_\
                  

$timestamp -> Prefix to name resulting files: $prefix

$timestamp -> Working directory: $wd
"

## START TIMER
start=$SECONDS

#### PIPELINE ####
create_wd && subsample && trim_adapters $wd$prefix"_sample.sorted.fastq" $wd$prefix"_sample.sorted.fastq" \
&& sort_file && reads_overlap && mt_reads_filt && first_assembly && gfa2fasta \
&& collecting_mt_reads $wd$prefix"_first_draft_asm.fasta" $input_file $wd$prefix"_align.bam" $wd$prefix"_collected_reads.fastq" \
&& trim_adapters $wd$prefix"_collected_reads.fastq" $wd$prefix"_collected_reads.fastq" && quality_control \
&& second_assembly && select_contig \
&& collecting_mt_reads $wd$prefix"_second_draft_asm.fasta" $input_file $wd$prefix"_align.bam" $wd$prefix"_collected_reads.fastq" \
&& trim_adapters $wd$prefix"_collected_reads.fastq" $wd$prefix"_collected_reads.fastq" \
&& quality_control && correct_reads && polish_asm

#correct_reads

echo ""
echo "### MITNANEX finished ###"
echo ""
echo "Final assembly is in" $wd$prefix"_final_draft_asm.fasta"
echo "Putative mitochondrial reads are in "$wd$prefix"_collected_reads.fastq"

## END TIMER
duration=$(( SECONDS - start ))
echo "" 
echo "$timestamp -> Elapsed time: $duration secs."
echo "" 