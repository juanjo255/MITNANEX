#!/bin/bash

## AESTHETICS
bold=$(tput bold)
normal=$(tput sgr0)
color_red='\033[0;91m'
color_cyan='\033[0;36m'
no_color='\033[0m'
timestamp=$(date -u +"%Y-%m-%d %T")

## IMPORT SCRIPTS
scripts_path="$exec_path/scripts"
source "$scripts_path/filter_by_quality.sh"
source "$scripts_path/map_reads.sh"
source "$scripts_path/variant_calling.sh"
#source "$scripts_path/haplogroup_class.sh"
source "$scripts_path/annotate_vcf.sh"
source "$scripts_path/annotate_mito.sh"
source "$scripts_path/mitomap.sh"
source "$scripts_path/plot_vcf.sh"


#DEFAULT
WD="."
output_folder="mitnanex_results"
threads="8"
minimap2_opts="-ax map-ont"
min_mapQ="30"
min_pruning=3
kmer_size="--kmer-size 15 --kmer-size 25"
medaka_model="r1041_e82_400bps_sup_variant_v5.0.0"
#haplogrep_trees="phylotree-rcrs@17.2"
#haplogrep_posible_trees=$("$exec_path/haplogrep3" trees)
top_hits="3"
keep_percent=50
min_length=500
max_length=2147483647
flye_preset="--nano-hq"
other_flye_opts=""
min_mean_quality=10
no_filter=1
ref_gff="$exec_path/refseqMT/chrMT_NC_012920.gff3"

## HELP MESSAGE
help() {
    echo "
    MITNANEX - MITochondrial NANopore reads EXtractor


    Usage: mitnanex.sh --reference genome.fasta --reads reads.fastq [options]

    ${bold}Options:${normal}
        -r, --reference    Reference organelle Genome [required].
        -g, --gff          Reference GFF3 file to liftover. [human.gff3].
        -i, --reads        Input file. [required].
        --ID               Use if your reference contains nuclear genome. The ID is the sentence inmidiately next to the '>' without spaces.
        -t, --threads      Threads. [$threads].
        --mm2              Minimap2 options.[$minimap2_opts].
        --WD               Working Directory to place results. [$WD]
        -o, --output       Directory name. [$output_folder].
        --mapq             Minimun mapping quality. [$min_mapQ].
        --min_length       Min read length. [$min_length].
        --max_length       Max read length. [reference length].
        --min_mean_quality Min average read quality. [$min_mean_quality].
        --no_filter        Do not filter reads by quality and length. [False].
        *                  Help.
    
    ${bold}GATK options:${normal}
    --max_assembly_region_size      Size of active region and assembly for variant discovery.[median read length].
    --min_pruning                   Min number of reads supporting edge during haplotype assembly.[3]. 
    -k, --kmer_size                     kmer size for building Debrujin graph for haplotype assembly. Comma separated values, e.g. value1,value2. [15,25].

    ${bold}Haplogrep options:${normal}
    --trees                PhyloTrees mt for Haplogrep3. comma separated. e.g. value1,value2. [$haplogrep_trees]. Possible options {$haplogrep_posible_trees}   
    --top_hits             Return the INT top hits. [$top_hits].
    
    ${bold}Flye options:${normal}
    --flye_preset            Sequencing technology. [$flye_preset].              
    --other_flye_opts        Other flye options. [" "].  
    "

    exit 1
}

## PARSE ARGUMENTS
ARGS=$(getopt -o "hr:i:g:t:k:o:" --long "help,reference:,gff:,reads:,output:,mm2:,WD:,threads:,ID:,min_pruning:,kmer_size:,
max_assembly_region_size:,trees:,top_hits:,keep_percent:,min_length:,flye_preset:,other_flye_opts:,min_mean_quality:,
max_length:,mapq:,no_filter" -n 'MITNANEX' -- "$@")
eval set -- "$ARGS"

while [ $# -gt 0 ];
do
    case "$1" in
    -r | --reference )
        ref_genome=$2
        shift 2
    ;;
     -g | --gff )
        ref_gff=$2
        shift 2
    ;;
    -i | --reads )
        reads=$2
        shift 2
    ;;
    -t | --threads )
        threads=$2
        shift 2
    ;;
    --ID )
        ID=$2
        shift 2
    ;;
    --mm2 )
        minimap2_opts=$2
        shift 2
    ;;
    --mapq )
        min_mapQ=$2
        shift 2
    ;;
    --WD )
        WD=$2
        shift 2
    ;;
    -o | --output )
        output_folder=$2
        shift 2
    ;;
    --max_assembly_region_size )
        median_read_len=$2
        shift 2
    ;;
    -k | --kmer_size )
        kmer_size=""
        IFS="," read -a kmers <<< "$2"
        for kmer in "${kmers[@]}";
        do
            kmer_size="$kmer_size --kmer-size $kmer"
        done
        shift 2
    ;;
    --min_pruning )
        min_pruning=$2
        shift 2
    ;;
    --trees )
        haplogrep_trees=$2
        shift 2
    ;;
    --top_hits )
        top_hits=$2
        shift 2
    ;;
    --min_length )
        min_length=$2
        shift 2
    ;;
    --max_length )
        max_length=$2
        shift 2
    ;;
    --min_mean_quality )
        min_mean_quality=$2
        shift 2
    ;;
    --flye_preset )
        flye_preset=$2
        shift 2
    ;;
    --other_flye_opts )
        other_flye_opts=$2
        shift 2
    ;;
    --no_filter )
        no_filter=0
        shift 1
    ;;
    -h | --help )
        help
    ;;
    -- ) 
        shift
        break
    ;;
    * )
        echo "${color_red}${color_red}ERROR${no_color}${no_color}: Invalid option $1. Use -h or --help to see options"
        exit 1
    ;;
    esac
done

# ARGUMENTS CHECK

if [ -z "$ref_genome" ];
then
  echo "[${color_red}ERROR${no_color}]: reference genome is required."
  echo " "
  help
fi

if [ -z "$reads" ];
then
  echo "[${color_red}ERROR${no_color}]: reads are required."
  echo " "
  help
fi

## Setting a correct working dir
if [ ${WD: -1} = / ];
then 
    WD=$WD$output_folder
else
    WD=$WD"/"$output_folder
fi

## Prefix name to use for the resulting files
if [ -z $prefix ];
then 
    prefix=$(basename $reads)
    prefix=${prefix%%.*}
fi

#FUNCTIONS WORKFLOW

create_wd(){
## CREATE WORKING DIRECTORY
    if ! [ -d $1 ]
    then
        mkdir -vp $1
    fi
}
custom_prints(){
    echo " "
    echo -e "$timestamp #########${color_cyan} $1 ${no_color}#########"
    echo " "
}

vars_definition(){
    # This is to keep in one place all outputs need it by other processes
    # In case you are starting a specific process you already have the variables defined
    
    # Map reads
    aln_file="$WD/$prefix.$ID.sorted.bam"
    MT_reads="$WD/$prefix.fastq"
    consensus_mitogenome="$WD/MT_genome.fasta"
    
    # Variant calling
    gatk_folder="$WD/VariantCall/gatk_mutect2"
    vcf_file="$gatk_folder/$prefix.$ID.gatk.filt.vcf"

    # annotate mito
    annot_folder="$WD/annotation/"
    annot_file="$annot_folder/$prefix.liftoff.gff"

    # Mitomap
    mitomap_out="$WD/variant_annot.mitomap.txt"
}


pipe_exec(){
    create_wd $WD

    ## Checking if there is only 1 reference genome
    if [ $(grep -c ">" $ref_genome) -gt 1 ];
        then
            echo "[${color_red}ERROR${no_color}]: Your reference genome contains more than 1 contig."
            exit 1
    else
        ID=$(grep -o "^>[^ ]*" $ref_genome | sed 's/>//g')
        max_length=$(seqkit fx2tab -l -n $ref_genome | cut -f 2)
    
        if [[ $? -eq 0 ]]; then
            if [[ $no_filter -eq 1 ]]; then
                filter_by_quality
            fi
            vars_definition
            map_reads
            #variant_calling
            #annotate_vcf
            #plot_vcf
            #annotate_mito
            #mitomap

        else
            echo "[${color_red}ERROR${no_color}] $?: Something wrong with reference genome."
            exit $?
        fi
    fi
}

pipe_exec
