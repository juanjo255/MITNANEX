import numpy as np
import time
from src.mitnanex import run
import pandas as pd

inicio = time.time()
a = run()
final = time.time()
tiempo = final - inicio
print("my code", tiempo)

inicio = time.time()
file = pd.read_csv("overlaps_talaro_18_07_2023.paf",delimiter='\t',header=None)
file [10] = file[10].astype(int)
file = file[file[10]>500]
file_grouped = (file.groupby(by=0))#.size()
print(file_grouped.size().sort_values(ascending=False)[:20].describe())
print(file_grouped.ngroups)
final = time.time()
tiempo = final - inicio
print("pandas", tiempo)
