from Bio import SeqIO

def convert_fq_to_fa (file:str, output:str) -> None:
    
    with open(file) as handle:
        count_sequences = SeqIO.convert(handle, "fastq", output, 'fasta')
        print(f"{count_sequences} sequences converted")