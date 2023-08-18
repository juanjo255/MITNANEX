import numpy as np

# file = pd.read_csv('overlaps_talaro_18_07_2023.paf',delimiter="\t", header=None)
# file['coverage'] = ((file[3]-file[2])/file[1])
# file['len'] = file[3]-file[2]
import time
def prueba ():
    
    a = []
    for i in range(1, 100000):
        a.append(i)
    inicio = time.time()
    [i for i in a if i%2==0]
    final = time.time()
    total = final-inicio
    print("lista",total)

    b = np.array([])
    for i in range(1, 100000):
        b = np.append(b,i)
    inicio = time.time()
    np.where(b%2==0)
    final = time.time()
    total = final-inicio
    print("array conditions",total)

prueba()