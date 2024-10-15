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
        -r, --reference    Reference Genome [required].
        -i, --reads        Input file. [required].
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
ARGS=$(getopt -o "hr:i:t:" --long "help,reference:,reads:,mm2:,threads:," -n 'MITNANEX' -- "$@")
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
        echo "ERROR: Invalid option"
    ;;
    esac
    break
done


if [ ${wd: -1} = / ];
then 
    wd=$wd$output_folder
else
    wd=$wd"/"$output_folder
fi

## PREFIX name to use for the resulting files
if [ -z $prefix ];
then 
    prefix=$(basename $reads)
    prefix=${prefix%%.*}
fi

map_reads(){
    
    minimap2 $minimap2_opts --split-prefix "temp" --secondary=no  $ref_genome $reads | samtools view -b --min-MQ $min_mapQ -F4 -T $ref_genome -o "$wd/$prefix.bam"
}

