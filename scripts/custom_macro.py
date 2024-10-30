from vcfstats.macros import continuous, categorical, aggregation

@categorical
def POS(variant):
    """Get the position for each variant as a numpy array."""
    return variant.POS    

