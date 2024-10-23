#!/bin/bash

bold=$(tput bold)
normal=$(tput sgr0)

#DEFAULT
WD="."
output_folder="mitnanex_results"
threads="8"
minimap2_opts="-ax map-ont"
min_mapQ="30"
min_pruning="3"
kmer_size="--kmer-size 15 --kmer-size 25"
medaka_model="r1041_e82_400bps_sup_variant_v5.0.0"
haplogrep_trees="phylotree-rcrs@17.2"
#haplogrep_posible_trees=$("$exec_path/haplogrep3" trees)
top_hits="3"
keep_percent=50
min_length=500
flye_preset="--nano-hq"
other_flye_opts=""
max_length=""
min_mean_quality=15

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
        --min_length       Min read length. [$min_length].
        --max_length       Max read length. [reference length].
        --min_mean_quality Min average read quality. [$min_mean_quality].
        *                  Help.
    
    ${bold}GATK options:${normal}
    --max_assembly_region_size      Size of active region and assembly for variant discovery.[median read length].
    --min_pruning                   Min number of reads supporting edge during haplotype assembly.[3]. 
    -k, --kmer_size                     kmer size for building Debrujin graph for haplotype assembly. Comma separated values, e.g. value1,value2. [15,25].

    ${bold}Haplogrep options:${normal}
    --trees                PhyloTrees mt for Haplogrep3. comma separated. e.g. value1,value2. [$haplogrep_trees]. Possible options {$haplogrep_posible_trees}   
    --top_hits             Return the INT top hits. [$top_hits].
    
    ${bold}Filtlong options:${normal}
    --keep_percent         Throw out the worst PERCENT(%) of reads. [$keep_percent].
    
    
    ${bold}Flye options:${normal}
    --flye_preset            Sequencing technology. [$flye_preset].              
    --other_flye_opts        Other flye options. [" "].  
    "

    exit 1
}

## PARSE ARGUMENTS
ARGS=$(getopt -o "hr:i:t:k:" --long "help,reference:,reads:,mm2:,WD:,threads:,ID:,min_pruning:,kmer_size:,
max_assembly_region_size:,trees:,top_hits:,keep_percent:,min_length:,flye_preset:,other_flye_opts:,min_mean_quality:,max_length:" -n 'MITNANEX' -- "$@")
eval set -- "$ARGS"

while [ $# -gt 0 ];
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
    --keep_percent )
        keep_percent=$2
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
    -h | --help )
        help
    ;;
    -- ) 
        shift
        break
    ;;
    * )
        echo "ERROR: Invalid option. Use -h or --help to see options"
        exit 1
    ;;
    esac
done

# ARGUMENTS CHECK

if [ -z "$ref_genome" ];
then
  echo "[ERROR]: reference genome is required."
  echo " "
  help
fi

if [ -z "$reads" ];
then
  echo "[ERROR]: reads are required."
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
    color_cyan='\033[0;36m'
    no_color='\033[0m'
    echo " "
    echo -e "#########${color_cyan} $1 ${no_color}#########"
    echo " "
}


map_reads(){
    
    custom_prints "Filter reads and mapping to reference"
    ## Seqkit output
    seqkit_output="$WD/$prefix.filtQ$min_mean_quality.fastq"
    seqkit seq -g --threads $threads --min-len $min_length --max-len $max_length $reads | \
    chopper -q $min_mean_quality -o $seqkit_output
    reads=$seqkit_output

    
    ## Map reads to reference
    ## GATK needs read groups. -R for that reason in Minimap2.
    
    ## Minimap2 output
    aln_file="$WD/$prefix.sorted.bam"
    minimap2  --split-prefix $prefix --secondary=no -g 1k $minimap2_opts $ref_genome $reads | \
    samtools view -@ $threads -b --min-MQ $min_mapQ -F2052 -T $ref_genome | \
    samtools sort -@ $threads -o $aln_file

    ## Assemble with flye to remove possible NUMTs 
    custom_prints "Assemble with MetaFlye to remove bad quality and some Numts "

    ## Output for first MT reads and flye
    MT_reads="$WD/$prefix_reads.$ID.fastq"
    flye_folder="$WD/flye_for_numts"    
    samtools fastq -@ $threads $aln_file -o $MT_reads
    flye -t $threads --meta $flye_preset $MT_reads -o $flye_folder 

    # Map to flye assembly
    minimap2 --secondary=no -R '@RG\tID:samplename\tSM:samplename' \ 
        $minimap2_opts $flye_folder"/assembly.fasta" $MT_reads | samtools view --threads $threads -b --min-MQ $min_mapQ -F2048 > $flye_folder"/aln_"$prefix".sorted.bam"
    
    # Retrieve the mitochondria in the flye assembly which is the one with the highest coverage. 
    contig_ID=$(sort -n -k3 $flye_folder"assembly_info.txt" | tail -n 1 | cut -f 1)
    samtools view -@ $threads -b -F2048 $flye_folder"/aln_"$prefix".sorted.bam" $contig_ID | samtools sort -@ $threads -o $aln_file
    samtools fastq -@ $threads $aln_file -o $MT_reads


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

    

    ## GATK output
    vcf_nofilt_file="$gatk_folder/$prefix.$ID.gatk.vcf"
    vcf_file="$gatk_folder/$prefix.$ID.gatk.filt.vcf"

    ## GATK
    gatk Mutect2 -R $ref_genome -L $ID --mitochondria-mode \
    --dont-use-soft-clipped-bases --max-assembly-region-size $median_read_len --min-pruning $min_pruning \
    $kmer_size -I $aln_file -O $vcf_nofilt_file && gatk FilterMutectCalls --mitochondria-mode -O $vcf_file \
    -R $ref_genome -V $vcf_nofilt_file

    
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

annotate_vcf(){
    #Annotate VFC with rCRS reference
    vcf_file_annotated="$gatk_folder/$prefix.$ID.gatk.annot.vcf"
    if [ -s $vcf_file ]; then
    bcftools annotate -a $vcf_file/HV.bed   $vcf_file -c "CHROM,FROM,TO,Hypervariable"  -h <(echo '##INFO=<ID=Hypervariable,Number=1,Type=String,Description="Hypervariable">') -o $vcf_file_annotated 
    bcftools annotate -a $vcf_file/HP.bed    $vcf_file -c "CHROM,FROM,TO,Homopolymer"  -h <(echo '##INFO=<ID=Homopolymer,Number=0,Type=Flag,Description="Homoloplymer">') -o $vcf_file_annotated 
    bcftools annotate -a $vcf_file/HS.bed    $vcf_file -c "CHROM,FROM,TO,Hotspot"  -h <(echo '##INFO=<ID=Hotspot,Number=0,Type=Flag,Description="Hotspot">')              -o $vcf_file_annotated 
    bcftools annotate -a $vcf_file/CDS.bed   $vcf_file -c "CHROM,FROM,TO,CDS" -h <(echo '##INFO=<ID=CDS,Number=1,Type=String,Description="CDS">')                         -o $vcf_file_annotated 
    bcftools annotate -a $vcf_file/RNR.bed   $vcf_file -c "CHROM,FROM,TO,RNR"  -h <(echo '##INFO=<ID=RNR,Number=1,Type=String,Description="rRNA">')                       -o $vcf_file_annotated 
    bcftools annotate -a $vcf_file/TRN.bed   $vcf_file -c "CHROM,FROM,TO,TRN"  -h <(echo '##INFO=<ID=TRN,Number=1,Type=String,Description="tRNA">')                       -o $vcf_file_annotated 
    bcftools annotate -a $vcf_file/DLOOP.bed $vcf_file -c "CHROM,FROM,TO,DLOOP"  -h <(echo '##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="DLOOP">')                    -o $vcf_file_annotated 
    fi
}

assemble_haplotype(){
    ## Assemble an individual mitogenome
    ## Maybe this could improved, so far the logic is to generate a consensus using a simple nucleotide frequency
    ## With filtlong we will keep the best reads, so consensus with flye is better set with the most frequent haplotype
    ## However, given that we are using very old people, idk how true this consensus will be.
    filtered_MT_reads="$WD/$prefix_reads.$ID.filtlong.fastq"
    filtlong --min_length $min_length --keep_percent $keep_percent $MT_reads > $filtered_MT_reads
    
    ## Define output
    flye_folder="$WD/assembly"
    flye -t $threads --meta --no-alt-contigs $other_flye_opts $flye_preset $filtered_MT_reads -o $flye_folder

}

pipe_exec(){
    create_wd $WD

    ## Checking if there is only 1 reference genome
    if [ $(grep -c ">" $ref_genome) -gt 1 ];
        then
            echo "[ERROR]: Your reference genome contains more than 1 contig."
            exit 1
    else
        ID=$(grep -o "^>[^ ]*" $ref_genome | sed 's/>//g')
        max_length=$(seqkit fx2tab -l -n $ref_genome | cut -f 2)
    fi

}
