#!/bin/bash

annotate_mito(){

    custom_prints "Lift genes from mitogenome reference to mitogenome of current sample"

    annot_folder="$WD/annotation/"
    create_wd $annot_folder
    annot_file="$annot_folder/$prefix.liftoff.gff"

    liftoff -g $ref_gff "$consensus_mitogenome" \
    "$exec_path/refseqMT/chrMT.fna" -o $annot_file

    echo "$timestamp [ATTENTION]: Mitochondrial consensus annotation file is at" $annot_file

}