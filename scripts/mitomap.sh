#!/bin/bash

mitomap() {
    ## Connect to MITOMAP API to annotate variants, haplogroup, clinic relevance
    custom_prints "Annotate variants with Mitomap"

    mitomap_out="$WD/variant_annot.mitomap.txt"

    # Convert VCF to TSV-like file
    #bcftools query -i 'FILTER="PASS"' -f '%CHROM\t%POS\t%REF\t%ALT' all_20240815_human.NC_012920.1.gatk.filt.anno.vcf > prueba.vcf

    $exec_path/scripts/mitomap.py $vcf_file $mitomap_out

    echo "$timestamp [ATTENTION]: Variant annotations are at" $mitomap_out

}