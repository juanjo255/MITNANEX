import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

## PAF FORMAT

# 0	string	Query sequence name
# 1	int	Query sequence length
# 2	int	Query start coordinate (0-based)
# 3	int	Query end coordinate (0-based)
# 4	char	‘+’ if query/target on the same strand; ‘-’ if opposite
# 5	string	Target sequence name
# 6	int	Target sequence length
# 7	int	Target start coordinate on the original strand
# 8	int	Target end coordinate on the original strand
# 9	int	Number of matching bases in the mapping
# 10	int	Number bases, including gaps, in the mapping
# 11	int	Mapping quality (0-255 with 255 for missing)


class cluster:
    """
    Store each cluster of reads, each cluster contains reads in the form of their id

    """

    def __init__(
        self, id_sequences: str, best_chain_score: int, longest_read: int
    ) -> None:
        self.id_sequences: list = list(id_sequences)
        self.best_chain_score: int = best_chain_score
        self.longest_read: int = longest_read
        self.total_bases: int = longest_read
        self.depth: int = self.total_bases // longest_read
        self.cluster_id = id(self)

    def add_sequence(self, id_sequence: str) -> None:
        """
        Add an ID read to the cluster

        Args:
            id_sequence (str): read ID
        """
        self.id_sequences.append(id_sequence)

    def compute_depth():
        pass

def cluster_reads():
    pass

def filter_alignment(alignment: str) -> cluster:
    """Take one line of a paf format given by minimap2, filter it and save it in a cluster
    Args:
        align (str): One line from paf file which contains the default paf format given
            by ava-ont function in minimap2.

    Returns:
        cluster: class where the read will be gathered otherwise it will initialized a new cluster.
    """
    ## One aligment information
    fields: list = alignment.split("\t")

    ## filter reads shorter than 500 pb
    length_threshold = 500
    if (int(fields[1].strip()) < length_threshold) or (
        int(fields[6]) < length_threshold
    ):
        return None

    ## filter only hits that where one read contains threshold
    ## NOTE: So far I do not take into account the difference between internal matches and containments
    ## Here containment is just the read aligned to other read at least threshold_contaiment
    threshold_containment = 0.8
    if fields[10] >= min(fields[1], fields[6]) * threshold_containment:
        cluster_reads()


# @profile # This to measure memory consumption
def open_paf(paf_file: str):
    file = open(paf_file, "r")
    alignment = file.readline()

    # while alignment:
    #    alignment = file.readline()


if __name__ == "__main__":
    chunksize = 1000
    # for chunk in pd.read_csv('overlaps_talaro_18_07_2023.paf', chunksize=chunksize, delimiter="\t"):
    #     print (chunk)
    #     break
    # file = np.genfromtxt('overlaps_talaro_18_07_2023.paf',delimiter='\t')
    open_paf("overlaps_talaro_18_07_2023.paf")
