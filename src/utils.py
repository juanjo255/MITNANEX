from Bio import SeqIO

def convert_fq_to_fa (fastq:str, output:str) -> None:
    with open(fastq) as handle:
        count_sequences = SeqIO.convert(handle, "fastq", output, 'fasta')
        print(f"{count_sequences} sequences converted to fasta")


def get_sequences_by_id(fasta:str, ids:list):
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in ids:
                #print(record.id)                
                yield (record.id, record.seq)

convert_fq_to_fa('/Users/jjpc/Documents/TESIS/MITNANEX_PROJECT/test/s_cervisae_CEN.PK113-7D_reads.sorted.fastq', '/Users/jjpc/Documents/TESIS/MITNANEX_PROJECT/test/s_cervisae_CEN.PK113-7D_reads.sorted.fasta')