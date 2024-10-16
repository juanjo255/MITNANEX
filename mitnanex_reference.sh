#!/bin/bash

#DEFAULT
wd="."
output_folder="mitnanex_results"
threads=4
minimap2_opts="-ax map-ont"
min_mapQ=30
min_pruning=3
kmer_size="--kmer-size 15 --kmer-size 25"

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
    
    GATK options:
    --max_assembly_region_size      Size of active region and assembly for variant discovery.[median read length].
    --min_pruning                   Min number of reads supporting edge during haplotype assembly.[3]. 
    -k, --kmer_size                     kmer size for building Debrujin graph for haplotype assembly. Comma separated values. [15,25].

    "
    exit 1
}

## PARSE ARGUMENTS
ARGS=$(getopt -o "hr:i:t:k:" --long "help,reference:,reads:,mm2:,threads:,ID:,min_pruning:,kmer_size:,max_assembly_region_size:," -n 'MITNANEX' -- "$@")
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
    --wd)
        wd=$2
        shift 2
    ;;
    -o | --output)
        output_folder=$2
        shift 2
    ;;
    --max_assembly_region_size)
        median_read_len=$2
        shift 2
    ;;
    -k | --kmer_size)
        kmer_size=""
        IFS="," read -a kmers <<< "$2"
        for kmer in "${kmers[@]}";
        do
            kmer_size="$kmer_size --kmer-size $kmer"
        done
        shift 2
    ;;
    --min_pruning)
        min_pruning=$2
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
if [ $(grep -c ">" $ref_genome) -gt "1" ];
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
    ## GATK needs read groups. -R for that reason.
    minimap2  --split-prefix "temp" --secondary=no -R '@RG\tID:samplename\tSM:samplename' $minimap2_opts $ref_genome $reads | \
    samtools view --threads $threads -b --min-MQ $min_mapQ -F4 -T $ref_genome | \
    samtools sort --threads $threads -o "$wd/$prefix.sorted.bam"
    aln_file="$wd/$prefix.sorted.bam"
}

select_contig (){
    ## Select organelle with to reference ID
    #samtools index "$wd/$prefix.bam" && samtools idxstats "$wd/$prefix.bam"
    samtools view -b "$wd/$prefix.sorted.bam" $ID > "$wd/$prefix.$ID.sorted.bam"
    aln_file="$wd/$prefix.$ID.sorted.bam"

    ## Separate mitogenome for future use
    seqkit grep -p $ID -o "$wd/$prefix.$ID.MT.fasta"
    ref_genome="$wd/refMT.$ID.fasta"
}

variant_calling() {
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
    samtools fastq $aln_file -o "$wd/$prefix_reads.$ID.fastq"
    MT_reads="$wd/$prefix_reads.$ID.fastq"

    ## GATK
    gatk Mutect2 -R $ref_genome -L $ID --mitochondria-mode \
    --dont-use-soft-clipped-bases --max-assembly-region-size $median_read_len --min-pruning $min_pruning \
    $kmer_size -I $aln_file -O "$wd/$prefix.$ID.vcf"

    ## Medaka
    medaka_variant -i $MT_reads -r $ref_genome -t $threads
    ## Sniffles


}

pipe_exec (){
    map_reads && echo " "
    if ! [ -z $ID ];then
        select_contig && echo " "
    else
        ID=$(grep -o "^>[^ ]*" $ref_genome | sed 's/>//g')
    fi

}