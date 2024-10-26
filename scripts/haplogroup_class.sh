#!/bin/bash

haplogroup_class(){
    
    custom_prints "Classify haplogroup"

    ## In case you are starting from here, you need this varibles
    gatk_folder="$WD/VariantCall/gatk_mutect2"
    vcf_file="$gatk_folder/$prefix.$ID.gatk.filt.vcf"
    
    haplogroup_folder="$WD/haplogroup"
    create_wd $haplogroup_folder

    ## Install trees
    "$exec_path/haplogrep3" install-tree $haplogrep_trees && echo " " || echo "${color_red} ERROR ${no_color} while downloading trees. Make sure .yaml has permissions and that you have internet"

    ## Classify 
    IFS="," read -a trees <<< "$haplogrep_trees"
        for tree in "${trees[@]}";
        do
            "$exec_path/haplogrep3" classify --tree=$tree --in $vcf_file --hits $top_hits \
                --extend-report --out "$haplogroup_folder/haplogrep3.$tree.txt" || echo "${color_red} ERROR ${no_color} Are the trees downloaded?"
        done

    ## SUMMARY RESULTS
    echo "$timestamp [ATTENTION]: The report with the top $top_hits closest haplogroups is at" "$haplogroup_folder/haplogrep3.$tree.txt"
    
}