
## Pairwise Sequence Aligment
class psa:
    
    def __init__(self, alignment_fields:list) -> None:
        self.query_id = alignment_fields[0]
        self.query_length = int(alignment_fields[1])
        self.begin_query = alignment_fields [2]
        self.end_query = alignment_fields[3]
        self.strand = alignment_fields[4]
        self.reference_id = alignment_fields[5]
        self.reference_length = int(alignment_fields[6])
        self.begin_reference = alignment_fields [7]
        self.end_reference = alignment_fields[8]
        self.match_bases = int(alignment_fields[9])
        self.align_length = int(alignment_fields[10])
        self.mapping_quality = int(alignment_fields[11])
        self.map_identity = self.match_bases / self.align_length
    