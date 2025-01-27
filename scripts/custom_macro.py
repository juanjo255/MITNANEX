from vcfstats.macros import continuous, categorical, aggregation
import numpy as np

@categorical
def POS(variant):
    """Get the position for each variant as a numpy array."""
    return variant.POS

    
@continuous
def AF(variant):
    """Get allele fractions of alternate alleles as a numpy array."""
    return np.concatenate([variant.format("AF").flatten()])

@continuous
def MBQ(variant):
    """Get median base quality by allele as a numpy array."""
    return [i for i in variant.INFO["MBQ"]]
