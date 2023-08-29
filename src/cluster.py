from .psa import psa
class cluster:
    """
    Store each cluster of reads, each cluster contains reads in the form of their id
    """

    def __init__(
        self, alignment:psa) -> None:
        self.id_sequences: list = [alignment.query_id, alignment.reference_id]
        # self.best_chain_score: int = alignment_fields[14]
        if alignment.reference_length > alignment.query_length:
            self.longest_read_length : int = alignment.reference_length
            self.longest_read_id:str = alignment.reference_id
        else:
            self.longest_read_length : int = alignment.query_length
            self.longest_read_id:str = alignment.query_id
        
        ## NOTE: This coverage is number of reads saved in this cluster
        self.coverage: int = 2 #self.total_bases // alignment_fields
        self.id_cluster = None

    def add_id(self, id_sequence: str) -> None:
        """
        Add an ID read to the cluster

        Args:
            id_sequence (str): read ID
        """
        self.id_sequences.append(id_sequence)
        self._update_coverage()

    def _update_coverage(self):
        """
        Update self.coverage 
        """
        self.coverage = len(self.id_sequences)

    def update_longest_read(self, read_length_longest:int, longest_read_id:str):
        """
        Update the length of the longest read in the cluster
        
        Args:
            read_length (int): length of the read to add
        """
        self.longest_read_length = read_length_longest
        self.longest_read_id = longest_read_id

    def set_cluster_id (self, id:int):
        self.id_cluster = id
        