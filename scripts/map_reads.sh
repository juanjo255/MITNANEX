#!/bin/bash

map_reads(){
    
    ## PRINT
    custom_prints "Mapping to reference"
    
    ## Map reads to reference
    ## GATK needs read groups. -R for that reason in Minimap2.
    
    ## Minimap2 output
    aln_file="$WD/$prefix.$ID.sorted.bam"
    minimap2 --secondary=no $minimap2_opts -g 1k  $ref_genome $reads | \
    samtools view -@ $threads -b --min-MQ $min_mapQ -F4,2048 -T $ref_genome | \
    samtools sort -@ $threads > $aln_file

    ## PRINT   
    ## Assemble with flye to remove possible NUMTs 
    custom_prints "Assemble with MetaFlye to remove bad quality and some Numts "

    ## Output for first MT reads and assemble with flye
    MT_reads="$WD/$prefix.fastq"
    
    flye_folder="$WD/flye_for_numts"    
    samtools fastq -@ $threads $aln_file > $MT_reads
    filtlong  --keep_percent $keep_percent --length_weight 10 $MT_reads > "$WD/$prefix.filtlong.fastq"
    MT_reads="$WD/$prefix.filtlong.fastq"
    rm $aln_file
    flye -t $threads --meta $flye_preset $MT_reads -o $flye_folder 

    # Map to flye assembly
    minimap2 --secondary=no $minimap2_opts -k 25 $flye_folder"/assembly.fasta" $MT_reads | \
    samtools view --threads $threads -b --min-MQ $min_mapQ -F4,2048 | samtools sort  -@ $threads > $flye_folder"/aln_"$prefix".sorted.bam"
    samtools index -@ $threads $flye_folder"/aln_"$prefix".sorted.bam"

    ## PRINT
    custom_prints "Retrieve mitochondria and remap reads"
    # Retrieve the mitochondria in the flye assembly which is the one with the highest depth or longest.
    # In assembly info field 3 is depth and 2 is lenght
    contig_ID=$(sort -n -k2 $flye_folder"/assembly_info.txt" | tail -n 1 | cut -f 1)

    ## Save mitogenome flye consensus
    consensus_mitogenome="$WD/MT_genome.fasta"
    seqkit grep -p $contig_ID "$flye_folder/assembly.fasta" > $consensus_mitogenome


    # To avoid confusion I remove reads used during assembly \
    # as I will keep only reads remapped to consensus
    rm $MT_reads
    
    ## Retrieve reads which mapped to the consensus_mitogenome 
    MT_reads="$WD/$prefix.fastq"
    samtools view  -@ $threads -b -F4 $flye_folder"/aln_"$prefix".sorted.bam" $contig_ID >  $flye_folder"/aln_"$prefix"_$contig_ID.bam"
    samtools sort  -@ $threads $flye_folder"/aln_"$prefix"_$contig_ID.bam" > $flye_folder"/aln_"$prefix"_$contig_ID.sorted.bam"
    samtools fastq -@ $threads $flye_folder"/aln_"$prefix"_$contig_ID.sorted.bam" > $MT_reads

    ## Removing unneeded files
    #rm $flye_folder"/aln_"$prefix"_$contig_ID.bam"
    #rm $flye_folder"/aln_"$prefix".sorted.bam"

    ## Final align file for variant calling 
    minimap2 --secondary=no -R '@RG\tID:samplename\tSM:samplename' $minimap2_opts $ref_genome $MT_reads | \
    samtools view -@ $threads -b -F4 -T $ref_genome | samtools sort -@ $threads > $aln_file
    samtools index -@ $threads $aln_file

    ## PRINT OUTPUT SUMMARY
    echo "$timestamp [ATTENTION]: Consensus mitogenome is at" $consensus_mitogenome
    echo "$timestamp [ATTENTION]: Mitochondrial reads are at: " $MT_reads
    echo "$timestamp [ATTENTION]: Mitochondrial reads mapped mitochondrial reference $(basename $ref_genome) at: " $aln_file
    echo "$timestamp [ATTENTION]: Mitochondrial reads mapped mitochondrial consensus at: " $flye_folder"/aln_"$prefix"_$contig_ID.sorted.bam"
    num_mapped_reads=$(samtools view -c $aln_file)
    echo "$timestamp Number of reads mapped: $aln_file"
}