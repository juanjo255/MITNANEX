import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

chunksize = 1000
for chunk in pd.read_csv('', chunksize=chunksize):
    print (chunk.head())
    break
