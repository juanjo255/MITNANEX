import pandas as pd
from Bio import SeqIO


def select_contig(flye_metadata: str, fasta: str = None) -> str:
    """Select the contig to use to use in minimap2 continue summoning reads for a better assembly.

    Returns:
        str: Sequence of the selected contig.
    """
    flye_metadata_df = pd.read_csv(flye_metadata, delimiter="\t", header=0)
    flye_metadata_df.sort_values(by=["length", "cov."], inplace=True)

    ## Check if there is a circular genome. If so and there are more than one, select the longest
    if flye_metadata_df.empty:
        print("Not assembly found!!")
        raise (ValueError("Not assembly found!"))
    elif not (flye_metadata_df[flye_metadata_df["circ."] == "Y"].empty):
        print("A circular genome found. Selecting... ")
        flye_metadata_df = flye_metadata_df[
            flye_metadata_df["circ."] == "Y"
        ].sort_values(by="length")
        return get_contig_sequence(flye_metadata_df["#seq_name"][0], fasta)
    else:
        return get_contig_sequence(flye_metadata_df["#seq_name"][0], fasta)


select_contig("/Users/jjpc/flye_s_cervisae_CEN_PK113-7D/assembly_info_2.txt")


def get_contig_sequence(contig: str, fasta: str) -> str:
    """look for the sequence of the selected contig

    Args:
        contig (str): selected contig.
        fasta (str): file returned by flye.

    Returns:
        str: Selected contig sequence.
    """
    with open(fasta) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id == contig:
                return record.seq
            else:
                raise (
                    ValueError(
                        "Something weird happend. I did not find the sequence of the contig."
                    )
                )
