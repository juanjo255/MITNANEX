from Bio import SeqIO

def convert_fq_to_fa (fastq:str, output:str) -> None:
    output = "".join(output.split(".")) + '.fasta'
    with open(fastq) as handle:
        count_sequences = SeqIO.convert(handle, "fastq", output, 'fasta')
        print(f"{count_sequences} sequences converted to fasta")
    return output

def get_sequences_by_id(fasta:str, ids:list) -> tuple:
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in ids:
                yield (record.id, record.seq)
    return
def write_fasta(fasta:str, sequences_ids:set, output:str):
    final_records = list()
    n=0
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if (record.id).strip() in sequences_ids:
                n+=1
                final_records.append(record)
                sequences_ids.remove(record.id)
    print(f"{n} reads retrieved")
    SeqIO.write(sequences=final_records, handle=output, format='fasta')