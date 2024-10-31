#!/bin/bash

map_reads(){
    
    ## PRINT
    custom_prints "Mapping to reference"
    
    ## Map reads to reference
    ## GATK needs read groups. -R for that reason in Minimap2.
    
    ## Minimap2 output
    aln_file="$WD/$prefix.$ID.sorted.bam"
    minimap2 --secondary=no $minimap2_opts -g 1k  $ref_genome $reads | \
    samtools view -@ $threads -b --min-MQ $min_mapQ -F2052 -T $ref_genome | \
    samtools sort -@ $threads > $aln_file

    ## PRINT   
    ## Assemble with flye to remove possible NUMTs 
    custom_prints "Assemble with MetaFlye to remove bad quality and some Numts "

    ## Output for first MT reads and assemble with flye
    MT_reads="$WD/$prefix.fastq"
    flye_folder="$WD/flye_for_numts"    
    samtools fastq -@ $threads $aln_file > $MT_reads
    rm $aln_file
    flye -t $threads --meta $flye_preset $MT_reads -o $flye_folder 

    # Map to flye assembly
    minimap2 --secondary=no $minimap2_opts -k 25 $flye_folder"/assembly.fasta" $MT_reads | \
    samtools view --threads $threads -b --min-MQ $min_mapQ -F2052 | samtools sort  -@ $threads > $flye_folder"/aln_"$prefix".sorted.bam"
    samtools index -@ $threads $flye_folder"/aln_"$prefix".sorted.bam"

    ## PRINT
    custom_prints "Retrieve mitochondria and remap reads"
    # Retrieve the mitochondria in the flye assembly which is the one with the highest coverage. 
    contig_ID=$(sort -n -k3 $flye_folder"/assembly_info.txt" | tail -n 1 | cut -f 1)

    ## Save mitogenome flye consensus
    consensus_mitogenome="$WD/MT_genome.fasta"
    seqkit grep -p $contig_ID "$flye_folder/assembly.fasta" > $consensus_mitogenome

    ## Retrieve reads which mapped to the consensus_mitogenome 
    samtools view  -@ $threads -b -F2052 $flye_folder"/aln_"$prefix".sorted.bam" $contig_ID >  $flye_folder"/aln_"$prefix"_$contig_ID.bam"
    samtools sort  -@ $threads $flye_folder"/aln_"$prefix"_$contig_ID.bam" > $flye_folder"/aln_"$prefix"_$contig_ID.sorted.bam"
    samtools fastq -@ $threads $flye_folder"/aln_"$prefix"_$contig_ID.sorted.bam" > $MT_reads

    ## Removing unneeded files
    #rm $flye_folder"/aln_"$prefix"_$contig_ID.bam"
    #rm $flye_folder"/aln_"$prefix".sorted.bam"

    ## Final align file for variant calling 
    minimap2 --secondary=no -R '@RG\tID:samplename\tSM:samplename' $minimap2_opts $ref_genome $MT_reads | \
    samtools view -@ $threads -b -F2052 -T $ref_genome | samtools sort -@ $threads > $aln_file
    samtools index -@ $threads $aln_file

    ## PRINT OUTPUT SUMMARY
    echo "$timestamp [ATTENTION]: Consensus mitogenome is at" $consensus_mitogenome
    echo "$timestamp [ATTENTION]: Mitochondrial reads are at: " $MT_reads
    echo "$timestamp [ATTENTION]: Mitochondrial reads mapped mitochondrial reference $(basename $ref_genome) at: " $aln_file
    echo "$timestamp [ATTENTION]: Mitochondrial reads mapped mitochondrial consensus at: " $flye_folder"/aln_"$prefix"_$contig_ID.sorted.bam"
    num_mapped_reads=$(samtools view -c $aln_file)
    echo "$timestamp Number of reads mapped: $aln_file"
}