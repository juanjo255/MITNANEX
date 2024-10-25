#!/bin/bash

#aesthetics
color='\033[0;91m'
no_color='\033[0m'

## Executable path
# FIXME: This part could generate problems, I am trying just to get the executable path to locate other folders
if which mitnanex.sh > /dev/null 2>&1; then
    ## This will not work for a while, I guess. 
    ## This is thought for a wet dream of using in anaconda or nextflow
    exec_path=$(grep -o ".*/" $(which mitnanex.sh))
else
    if [ "$(grep -o "/" <<< $0 | wc -l)" -gt 0 ]; then
        exec_path=$(grep -o ".*/" <<< $0)
    else
        exec_path="./"
    fi
    export exec_path=$exec_path
fi

## Help message
help() {
    echo -e "
    MITNANEX - MITochondrial NANopore reads EXtractor

    Version: 0.1

    Usage: mitnanex.sh [subcommand][options] FASTQ

    subcommands:
        ${color} denovo ${no_color}  - (On Development) Run denovo strategy using kmer analysis.
        ${color} reference ${no_color} - Run using a reference genome.
    
    
    Authors:
    Juan Jose Picon Cossio
    Javier Correa Alvarez
    Gustavo Gamez de las Armas
    "
    exit 1
}
while true;do
    case $1 in
    denovo)
        shift 1
        bash "$exec_path/mitnanex_denovo.sh" $@
    ;;
    reference)
        shift 1
        bash "$exec_path/mitnanex_reference.sh" $@
    ;;
    --help | -h)
        help 
    ;;
    *)
        echo "ERROR: Invalid option"
        help
    ;;
    esac
    break
done
