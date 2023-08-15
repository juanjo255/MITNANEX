import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class cluster ():
    """
    Store each cluster of reads, each cluster contains reads in the form of their id
    
    """
    def __init__(self, id_sequences:str, best_chain_score:int, longest_read:int) -> None:
        self.id_sequences:list = list(id_sequences)
        self.best_chain_score:int = best_chain_score
        self.longest_read:int = longest_read
        self.total_bases:int = longest_read
        self.depth:int = self.total_bases//longest_read

    def add_sequence (self, id_sequence:str) -> None:
        """
        Add an ID read to the cluster

        Args:
            id_sequence (str): read ID
        """
        self.id_sequences.append(id_sequence)

    def compute_depth ():
        pass

def get_alignment (align:str) -> cluster:
    """ Get one line of a paf format given by minimap2

    Args:
        align (str): One line from paf file which contains the default paf format given by ava-ont function in minimap2  

    Returns:
        cluster: class where the read will be gathered otherwise it will initialized a new cluster
    """
    pass

#@profile # This to measure memory consumption
def open_paf (paf_file:str):
    file = open(paf_file, 'r')
    alignment = file.readline()
    
    #while alignment:
    #    alignment = file.readline()
    

if __name__ == '__main__':
    chunksize = 1000
    # for chunk in pd.read_csv('overlaps_talaro_18_07_2023.paf', chunksize=chunksize, delimiter="\t"):
    #     print (chunk)
    #     break
    # file = np.genfromtxt('overlaps_talaro_18_07_2023.paf',delimiter='\t')
    open_paf('overlaps_talaro_18_07_2023.paf')

