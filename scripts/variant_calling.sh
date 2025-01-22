#!/bin/bash

variant_calling() {

    custom_prints "Starting Variant calling"
    
    ## Create dirs
    gatk_folder="$WD/VariantCall/gatk_mutect2"
    create_wd $gatk_folder

    ## Variant calling with GATK
    
    if [ -z $median_read_len ];then
        median_read_len=$(cramino -t $threads $aln_file | grep "Median length" | cut -f 2)
        median_read_len=${median_read_len%%.*}
        median_read_len=$(($median_read_len + 1))
    fi
    
    #preprocessing files for tools

    ## Create index and dict
    #samtools faidx $ref_genome
    #gatk CreateSequenceDictionary -R $ref_genome

    

    ## GATK output
    vcf_nofilt_file="$gatk_folder/$prefix.$ID.gatk.vcf"
    vcf_file="$gatk_folder/$prefix.$ID.gatk.filt.vcf"

    ## GATK
    #--max-assembly-region-size $median_read_len
    gatk Mutect2 -R $ref_genome -L $ID --mitochondria-mode \
    --annotation "DepthPerAlleleBySample" --min-pruning $min_pruning \
    $kmer_size -I $aln_file -O $vcf_nofilt_file && \
    gatk FilterMutectCalls --mitochondria-mode -O $vcf_file -R $ref_genome -V $vcf_nofilt_file

    echo "$timestamp [ATTENTION]: The variant calling is at" $vcf_file

    
}