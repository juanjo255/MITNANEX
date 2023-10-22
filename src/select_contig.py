import pandas as pd
from Bio import SeqIO
import sys

def select_contig(flye_metadata: str, flye_asm_file: str) -> str:
    """Select the contig to use to use in minimap2 continue summoning reads for a better assembly.

    Returns:
        str: Sequence of the selected contig.
    """
    flye_metadata_df = pd.read_csv(flye_metadata, delimiter="\t", header=0)
    flye_metadata_df.sort_values(by=["cov.", "length"], inplace=True)

    ## Check if there is a circular genome. If so and there are more than one, select the longest
    if flye_metadata_df.empty:
        print("Not assembly found!!")
        raise (ValueError("Not assembly found!"))
    
    elif not (flye_metadata_df[flye_metadata_df["circ."] == "Y"].empty):
        print("A circular genome found!!")
        flye_metadata_df = flye_metadata_df[
            flye_metadata_df["circ."] == "Y"
        ].sort_values(by=["cov.", "length"])
        return get_contig_sequence(flye_metadata_df["#seq_name"][0], flye_asm_file)
    
    else:
        print("No circular genome found")
        return get_contig_sequence(flye_metadata_df["#seq_name"][0], flye_asm_file)


def get_contig_sequence(contig: str, flye_asm_file: str) -> str:
    """look for the sequence of the selected contig

    Args:
        contig (str): selected contig.
        flye_asm_file (str): file returned by flye.

    Returns:
        str: Selected contig sequence.
    """
    with open(flye_asm_file) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id == contig:
                print("Contig found: ",record.id)
                return record
            else:
                raise (
                    ValueError(
                        "Something weird happend. I did not find the sequence of the contig."
                    )
                )

if __name__ == "__main__":
    
    args = sys.argv
    try:
        path_flye_fasta = args[1] + "assembly_info.txt"
        assembly_file = args[1] + "assembly.flye_asm_file"
        output_file = args[2]
    except:
        raise(ValueError("The path to the flye results is missing!"))

    SeqIO.write(sequences=select_contig(path_flye_fasta, assembly_file), handle=output_file, format='fasta')
    