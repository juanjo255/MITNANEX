#!/bin/bash

mitomap() {
    ## Connect to MITOMAP API to annotate variants, haplogroup, clinic relevance
    custom_prints "Annotate variants with Mitomap"

    mitomap_out="$WD/variant_annot.mitomap.txt"

    # Convert VCF to TSV-like file
    bcftools query -i 'FILTER="PASS"' -f '%CHROM\t%POS\t%REF\t%ALT' $vcf_file > "$vcf_file.tsv"
    trans_vcf="$vcf_file.tsv"
    $exec_path/scripts/mitomap.py $trans_vcf $mitomap_out

    echo "$timestamp [ATTENTION]: Variant annotations are at" $mitomap_out

}