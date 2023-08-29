
def estimate_hash_table_size(path:str):
    """number of reads in fast file

    Args:
        path (str): location of the file
    """
    file = open(path, "r")
    line = file.readline().strip()
    size = 0
    while line:
        if line.startswith(('A','T','G','C')):
            size+=1
        # New alignment
        line = file.readline().strip()
    file.close()
    print(size)

