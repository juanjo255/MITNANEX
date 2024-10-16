#!/bin/bash

#DEFAULT
wd="."
output_folder="mitnanex_results"
threads=4
min_qual=0
minimap2_opts="-ax map-ont"
min_mapQ=30


## HELP MESSAGE
help() {
    echo "
    MITNANEX - MITochondrial NANopore reads EXtractor


    Usage: mitnanex.sh --reference genome.fasta --reads reads.fastq [options]

    Options:
        -r, --reference    Reference organelle Genome [required].
        -i, --reads        Input file. [required].
        --ID               Use if your reference contains nuclear genome. The ID is the sentence inmidiately next to the '>' without spaces.
        -t, --threads      Threads. [$threads].
        --mm2              Minimap2 options.[$minimap2_opts].
        --wd               Working Directory to place results. [$wd]
        -o, --output       Outout directory. [$output_folder].
        --mapq             Minimun mapping quality. [$min_mapQ].
        *                  Help.
    "
    exit 1
}

## PARSE ARGUMENTS
ARGS=$(getopt -o "hr:i:t:" --long "help,reference:,reads:,mm2:,threads:,ID," -n 'MITNANEX' -- "$@")
eval set -- "$ARGS"

while true;do
    case $1 in
    -r | --reference)
        ref_genome=$2
        shift 2
    ;;
    -i | --reads)
        reads=$2
        shift 2
    ;;
    --ID)
        ID=$2
        shift 2
    ;;
    --mm2)
        minimap2_opts=$2
        shift 2
    ;;
    --mapq)
        min_mapQ=$2
        shift 2
    ;;
    --help | -h)
        help 
    ;;
    --) 
        shift 1
    ;;
    *)
        echo "ERROR: Invalid option. Use -h or --help to see options"
    ;;
    esac
    break
done

# ARGUMENTS CHECK

## Checking if there is only 1 reference genome
if [ $(cat  $ref_genome | grep -c ">") -gt "1" ];
    then
        echo "[ERROR] Your reference genome contains more than 1 contig. Set --ID"
        exit 1
fi 

## Setting a correct working dir
if [ ${wd: -1} = / ];
then 
    wd=$wd$output_folder
else
    wd=$wd"/"$output_folder
fi

## Prefix name to use for the resulting files
if [ -z $prefix ];
then 
    prefix=$(basename $reads)
    prefix=${prefix%%.*}
fi



#FUNCTIONS WORKFLOW
map_reads(){
    ## Map reads to reference
    minimap2 $minimap2_opts --split-prefix "temp" --secondary=no  $ref_genome $reads | \
    samtools view --threads $threads -b --min-MQ $min_mapQ -F4 -T $ref_genome | \
    samtools sort --threads $threads -o "$wd/$prefix.sorted.bam"
    aln_file="$wd/$prefix.sorted.bam"
}

select_contig (){
    ## Select organelle according to reference ID
    #samtools index "$wd/$prefix.bam" && samtools idxstats "$wd/$prefix.bam"
    samtools view -b "$wd/$prefix.sorted.bam" $ID > "$wd/$prefix.$ID.sorted.bam"
    aln_file="$wd/$prefix.$ID.sorted.bam"
}

pipe_exec (){
    map_reads && echo " "
    if ! [ -z $ID ];then
        select_contig && echo " "
    fi
    
}