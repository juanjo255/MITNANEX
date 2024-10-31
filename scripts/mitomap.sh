#!/bin/bash

mitomap() {
    ## Connect to MITOMAP API to annotate variants, haplogroup, clinic relevance
    custom_prints "Annotate variants with Mitomap"

    mitomap_out="$WD/variant_annot.mitomap.txt"
    $exec_path/scripts/mitomap.py $vcf_file $mitomap_out

    echo "$timestamp [ATTENTION]: Variant annotations are at" $mitomap_out

}