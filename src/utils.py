from Bio import SeqIO

def convert_fq_to_fa (fastq:str, output:str) -> None:
    with open(fastq) as handle:
        count_sequences = SeqIO.convert(handle, "fastq", output, 'fasta')
        print(f"{count_sequences} sequences converted to fasta")
    return

def get_sequences_by_id(fasta:str, ids:list) -> tuple:
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in ids:
                yield (record.id, record.seq)
    return
def write_fasta(fasta:str, sequences_ids:list, output:str):
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in sequences_ids:
                SeqIO.write(sequences=record, handle=output, format='fasta')
    return