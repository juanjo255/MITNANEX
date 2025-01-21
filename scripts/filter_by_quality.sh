#!/bin/bash

filter_by_quality(){
    
    ## PRINT
    custom_prints "Filter reads by: Qual:$min_mean_quality MinLen: $min_length MaxLen: $max_length"
    
    ## Seqkit output
    chopper_output="$WD/$prefix.filtQ$min_mean_quality.fastq"
    chopper --threads $threads -q $min_mean_quality --minlength $min_length --maxlength $max_length < $reads > $chopper_output
    reads=$chopper_output
}