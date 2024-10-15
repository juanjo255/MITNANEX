#!/bin/bash

#DEFAULT
wd="."
output_folder="mitnanex_results"
threads=4
min_qual=0


## Help message
mitnanex_help() {
    echo "
    MITNANEX - MITochondrial NANopore reads EXtractor

    Author:
    Juan Picon Cossio

    Version: 0.1

    Usage: mitnanex.sh [options] FASTQ

    Options:
        -i        Input file. [required].
        -t        Threads. [4].
        -m        Min-len. Filter reads by minimun length. Read seqkit seq documentation. [-1].
        -M        Max-len. Filter reads by maximun length. Read seqkit seq documentation. [-1].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [$wd].
        -s        Mapping identity. Minimun identity between two reads to be store in the same cluster.[0.6].
        -q        Min mapping quality (>=). This is for samtools. [0].
        -f        Flye mode. [--nano-hq]
        -k        Percentage of reads to keep during filter with filtlong. [80]. 
        *         Help.
    "
    exit 1
}


map_reads(){
    # Input
    ## $1=reference
    ## $2=query
    ## $3=output
    
    minimap2 -ax map-ont --split-prefix "temp_prefix" --secondary=no  $1 $2 | samtools view -b --min-MQ $min_qual -F4 -T $1 -o $3

}
