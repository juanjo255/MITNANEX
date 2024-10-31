#!/bin/bash

plot_vcf(){
    ## PLOT VCF STATS

    custom_prints "Plot VCF statistics"

    vcfplotsdir="$WD/vcf_plots"
    create_wd $vcfplotsdir

    
    vcfstats --vcf $vcf_file --outdir $vcfplotsdir \
    --config "$exec_path/scripts/plots_vcf.toml" \
    --macro "$exec_path/scripts/custom_macro.py"

    echo "$timestamp [ATTENTION]: VCF plots are at" $vcfplotsdir

}