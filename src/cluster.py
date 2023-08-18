class cluster:
    """
    Store each cluster of reads, each cluster contains reads in the form of their id

    """

    def __init__(
        self, alignment_fields) -> None:
        self.id_sequences: list = list(alignment_fields[0],alignment_fields[5])
        # self.best_chain_score: int = alignment_fields[14]
        # self.longest_read: int = alignment_fields
        # self.total_bases: int = alignment_fields
        # self.coverage: int = self.total_bases // alignment_fields

    def add_id(self, id_sequence: str) -> None:
        """
        Add an ID read to the cluster

        Args:
            id_sequence (str): read ID
        """
        self.id_sequences.append(id_sequence)