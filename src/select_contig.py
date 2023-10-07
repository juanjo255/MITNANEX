from Bio import SeqIO

def select_contig (flye_metadata:str, fasta:str):
    """ Sometimes we can get more than one contig because of the 

    Args:
        fasta (str): _description_
        size_overlap (int, optional): _description_. Defaults to 1000.
        ksize (int, optional): _description_. Defaults to 3.

    Returns:
        _type_: _description_
    """
    with open(flye_metadata) as flye_metadata_handle:
        with open(fasta) as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                pass

            return 