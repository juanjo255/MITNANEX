from Bio import SeqIO
from sourmash import MinHash


def check_circularization(fasta: str, size_overlap: int = 1000, ksize=3):
    """To check if the genome is circular I collect seeds between the start and the final part of the
    genome and compared the jaccard similarity

    Args:
        fasta (str): _description_
        size_overlap (int, optional): _description_. Defaults to 1000.
        ksize (int, optional): _description_. Defaults to 3.

    Returns:
        _type_: _description_
    """
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_size = len(record.seq)
            start = record.seq[0:size_overlap]
            end = record.seq[seq_size + 1 - size_overlap : seq_size]
            start_signature = MinHash(n=0, ksize=ksize, scaled=1, track_abundance=True)
            start_signature.add_sequence(str(start))
            end_signature = MinHash(n=0, ksize=ksize, scaled=1, track_abundance=True)
            end_signature.add_sequence(str(end))

            return start_signature, end_signature


start_seq, end_seq = check_circularization(
    fasta="/Users/jjpc/flye_s_cervisae_CEN_PK113-7D/assembly.fasta",
    size_overlap=1000,
    ksize=3,
)
## check similary of the starting and final sequence
jaccard_simil = end_seq.jaccard(start_seq)
if jaccard_simil >= 0.9:
    print("Probably, we have circularized sequence!")
else:
    print("Jummm, genome seems incomplete! :( ")
