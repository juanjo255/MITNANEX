#!/bin/bash
annotate_vcf(){

    custom_prints "Annotate variants with reference features context"

    ## Temp file to annotate. The final annotation file will be $vcf_file
    vcf_file_annot="$gatk_folder/$prefix.$ID.gatk.filt.temp.vcf"

    ## Keep old name
    vcf_file_old=$vcf_file

    #Annotate VFC with rCRS reference
    reference_annot=$exec_path"/refseqMT"

    bcftools annotate -a "$reference_annot/HV.bed"   $vcf_file -c "CHROM,FROM,TO,Hypervariable"  -h <(echo '##INFO=<ID=Hypervariable,Number=1,Type=String,Description="Hypervariable">') > $vcf_file_annot
    cp $vcf_file_annot $vcf_file 
    bcftools annotate -a "$reference_annot/HP.bed"    $vcf_file -c "CHROM,FROM,TO,Homopolymer"  -h <(echo '##INFO=<ID=Homopolymer,Number=0,Type=Flag,Description="Homoloplymer">') > $vcf_file_annot 
    cp $vcf_file_annot $vcf_file
    bcftools annotate -a "$reference_annot/HS.bed"    $vcf_file -c "CHROM,FROM,TO,Hotspot"  -h <(echo '##INFO=<ID=Hotspot,Number=0,Type=Flag,Description="Hotspot">')              > $vcf_file_annot 
    cp $vcf_file_annot $vcf_file
    bcftools annotate -a "$reference_annot/CDS.bed"   $vcf_file -c "CHROM,FROM,TO,CDS" -h <(echo '##INFO=<ID=CDS,Number=1,Type=String,Description="CDS">')                         > $vcf_file_annot 
    cp $vcf_file_annot $vcf_file
    bcftools annotate -a "$reference_annot/RNR.bed"   $vcf_file -c "CHROM,FROM,TO,RNR"  -h <(echo '##INFO=<ID=RNR,Number=1,Type=String,Description="rRNA">')                       > $vcf_file_annot 
    cp $vcf_file_annot $vcf_file
    bcftools annotate -a "$reference_annot/TRN.bed"   $vcf_file -c "CHROM,FROM,TO,TRN"  -h <(echo '##INFO=<ID=TRN,Number=1,Type=String,Description="tRNA">')                       > $vcf_file_annot 
    cp $vcf_file_annot $vcf_file
    bcftools annotate -a "$reference_annot/DLOOP.bed" $vcf_file -c "CHROM,FROM,TO,DLOOP"  -h <(echo '##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="DLOOP">')                    > $vcf_file_annot 
    
    ## Change name to old name
    mv $vcf_file_annot $vcf_file_old

    echo "$timestamp [ATTENTION]: The annotated VCF is at" $vcf_file
}