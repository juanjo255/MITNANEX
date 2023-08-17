from typing import Any
import pandas as pd

# file = pd.read_csv('overlaps_talaro_18_07_2023.paf',delimiter="\t", header=None)
# file['coverage'] = ((file[3]-file[2])/file[1])
# file['len'] = file[3]-file[2]
import time

init1=time.time()
a = [i for i in range(1,10000000)]
if 300 in a:
    resultado = True
final1=time.time()
time1 = final1-init1
print(time1)

init2=time.time()
b = set(range(1,10000000))
if 300 in b:
    result = True
final2=time.time()
time2 = final2-init2
print(time2)
