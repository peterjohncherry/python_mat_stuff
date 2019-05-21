import numpy as np

def text_to_float(filename):
    with open(filename, 'r') as mat_file:
        input_matrix = [[float(num) for num in line.split(' ')] for line in mat_file]
    return np.array(input_matrix)
