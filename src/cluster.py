from psa import psa
class cluster:
    """
    Store each cluster of reads, each cluster contains reads in the form of their id
    """

    def __init__(
        self, alignment:psa) -> None:
        self.id_sequences: list = [alignment.query_id, alignment.reference_id]
        # self.best_chain_score: int = alignment_fields[14]
        # self.longest_read: int = alignment_fields
        # self.total_bases: int = alignment_fields
        self.coverage: int = 2 #self.total_bases // alignment_fields

    def add_id(self, id_sequence: str) -> None:
        """
        Add an ID read to the cluster

        Args:
            id_sequence (str): read ID
        """
        self.id_sequences.append(id_sequence)
        self._update_depth()

    def _update_depth(self):
        """
        Update self.coverage 
        """
        self.coverage = len(self.id_sequences)