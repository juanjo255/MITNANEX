## Parse mitomaster output, generate some plots
import pandas as pd

if __name__ == "__main__":
    file=pd.read_csv("/Users/jjpc/PiconCossio/UNI-SEQ/COL-HUMAN-PROJECT/20240918_human/variant_annot.mitomap.txt",
                    delimiter="\t")
    