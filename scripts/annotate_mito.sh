#!/bin/bash

annotate_mito(){

    custom_prints "Lift genes from mitogenome reference to mitogenome of current sample"

    annotation_file="$WD/annotation/$prefix.liftoff.gff"

    liftoff -g "$exec_path/refseqMT/chrMT_NC_012920.gff3" "$consensus_mitogenome" \
    "$exec_path/refseqMT/chrMT.fna" -o $annotation_file

    echo "$timestamp [ATTENTION]: Mitochondrial consensus annotation file is at" $annotation_file

}