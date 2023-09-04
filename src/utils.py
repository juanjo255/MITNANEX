from Bio import SeqIO

def convert_fq_to_fa (fastq:str, output:str) -> None:
    
    with open(fastq) as handle:
        count_sequences = SeqIO.convert(handle, "fastq", output, 'fasta')
        print(f"{count_sequences} sequences converted")


def get_sequences_by_id(fasta:str, ids:list):
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in ids:
                #print(record.id)                
                yield (record.id, record.seq)