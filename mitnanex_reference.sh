#!/bin/bash

bold=$(tput bold)
normal=$(tput sgr0)

#DEFAULT
WD="."
output_folder="mitnanex_results"
threads="4"
minimap2_opts="-ax map-ont"
min_mapQ="20"
min_pruning="3"
kmer_size="--kmer-size 15 --kmer-size 25"
medaka_model="r1041_e82_400bps_sup_variant_v5.0.0"
haplogrep_trees="phylotree-rcrs@17.2"
haplogrep_posible_trees=$("$exec_path/haplogrep3" trees)
top_hits="3"

## HELP MESSAGE
help() {
    echo "
    MITNANEX - MITochondrial NANopore reads EXtractor


    Usage: mitnanex.sh --reference genome.fasta --reads reads.fastq [options]

    ${bold}Options:${normal}
        -r, --reference    Reference organelle Genome [required].
        -i, --reads        Input file. [required].
        --ID               Use if your reference contains nuclear genome. The ID is the sentence inmidiately next to the '>' without spaces.
        -t, --threads      Threads. [$threads].
        --mm2              Minimap2 options.[$minimap2_opts].
        --WD               Working Directory to place results. [$WD]
        -o, --output       Outout directory. [$output_folder].
        --mapq             Minimun mapping quality. [$min_mapQ].
        *                  Help.
    
    ${bold}GATK options:${normal}
    --max_assembly_region_size      Size of active region and assembly for variant discovery.[median read length].
    --min_pruning                   Min number of reads supporting edge during haplotype assembly.[3]. 
    -k, --kmer_size                     kmer size for building Debrujin graph for haplotype assembly. Comma separated values, e.g. value1,value2. [15,25].

    ${bold}Medaka options:${normal}
    -m, model              Medaka model. [$medaka_model]

    ${bold}Haplogrep options:${normal}
    --trees                PhyloTrees mt for Haplogrep3. comma separated. e.g. value1,value2. [$haplogrep_trees]. Possible options {$haplogrep_posible_trees}   
    --top_hits             Return the INT top hits. [$top_hits].
    "
    exit 1
}

## PARSE ARGUMENTS
ARGS=$(getopt -o "hr:i:t:k:" --long "help,reference:,reads:,mm2:,WD:,threads:,ID:,min_pruning:,kmer_size:,max_assembly_region_size:,trees:,top_hits:," -n 'MITNANEX' -- "$@")
eval set -- "$ARGS"
while true;
do
    case "$1" in
    -r | --reference )
        ref_genome=$2
        shift 2
    ;;
    -i | --reads )
        reads=$2
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
    -h | --help )
        help
    ;;
    -- ) 
        shift
        break
    ;;
    * )
        echo "ERROR: Invalid option. Use -h or --help to see options"
        break
    ;;
    esac
done

# ARGUMENTS CHECK

## Checking if there is only 1 reference genome
if [ $(grep -c ">" $ref_genome) -gt "1" ];
    then
        echo "[ERROR] Your reference genome contains more than 1 contig. Set --ID"
        exit 1
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

map_reads(){
    ## Map reads to reference
    ## GATK needs read groups. -R for that reason.
    minimap2  --split-prefix "temp" --secondary=no -R '@RG\tID:samplename\tSM:samplename' $minimap2_opts $ref_genome $reads | \
    samtools view --threads $threads -b --min-MQ $min_mapQ -F4 -T $ref_genome | \
    samtools sort --threads $threads -o "$WD/$prefix.sorted.bam"
    aln_file="$WD/$prefix.sorted.bam"
}

select_contig(){
    ## Select organelle with to reference ID
    #samtools index "$WD/$prefix.bam" && samtools idxstats "$WD/$prefix.bam"
    samtools view -b "$WD/$prefix.sorted.bam" $ID > "$WD/$prefix.$ID.sorted.bam"
    aln_file="$WD/$prefix.$ID.sorted.bam"

    ## Separate mitogenome for future use
    seqkit grep -p $ID -o "$WD/$prefix.$ID.MT.fasta"
    ref_genome="$WD/refMT.$ID.fasta"
}

variant_calling() {

    ## Create dirs
    var_call_folder="$WD/VariantCall/"
    gatk_folder="$WD/$var_call_folder/gatk_mutect2"
    medaka_folder="$WD/$var_call_folder/medaka"

    create_wd $var_call_folder
    create_wd $gatk_folder
    create_wd $medaka_folder

    ## Variant calling with GATK and Medaka
    if [ -z $median_read_len ];then
        median_read_len=$(cramino $aln_file | grep "Median length" | cut -f 2)
    fi
    
    #preprocessing files for tools

    ## Create index and dict
    samtools index $aln_file
    samtools faidx $ref_genome
    gatk CreateSequenceDictionary -R $ref_genome

    ## BAM to fastq
    samtools fastq $aln_file -o "$WD/$prefix_reads.$ID.fastq"
    MT_reads="$WD/$prefix_reads.$ID.fastq"

    ## GATK
    gatk Mutect2 -R $ref_genome -L $ID --mitochondria-mode \
    --dont-use-soft-clipped-bases --max-assembly-region-size $median_read_len --min-pruning $min_pruning \
    $kmer_size -I $aln_file -O "$gatk_folder/$prefix.$ID.gatk.vcf" && gatk FilterMutectCalls --mitochondria-mode -O "$gatk_folder/$prefix.$ID.gatk.filt.vcf" \
    -R $ref_genome -V "$gatk_folder/$prefix.$ID.gatk.vcf"

    ## Medaka
    medaka_variant -t $threads -m $medaka_model -i $MT_reads -r $ref_genome  -o $medaka_folder

    ## RETURN
    vcf_file="$gatk_folder/$prefix.$ID.gatk.vcf"
}

haplogroup_class(){
    haplogroup_folder="$WD/haplogroup/"
    create_wd $haplogroup_folder

    ## Install trees
    "$exec_path/haplogrep3" install-tree $haplogrep_trees && echo " " || echo "Error while downloading trees. Make sure .yaml has permissions" exit 1

    ## Classify 
    IFS="," read -a trees <<< "$haplogrep_trees"
        for tree in "${trees[@]}";
        do
            "$exec_path/haplogrep3" classify --tree=$tree --in $vcf_file --hits $top_hits \
                --extend-report --out "$haplogroup_folder/haplogrep3.$tree"
        done

}

pipe_exec(){
    create_wd $WD
    map_reads && echo " "
    if ! [ -z $ID ];then
        select_contig && echo " "
    else
        ID=$(grep -o "^>[^ ]*" $ref_genome | sed 's/>//g')
    fi

}

