import numpy as np
import cmath
import davidson as dv
import mat_utils as util
import plot_tools as pt
import matplotlib.pyplot as plt
import read_mat

def main():
    # define dimensions of the array
    nrows = 6
    ncols = 6
    x_mat = util.generate_random_matrix(nrows, ncols)  # type: List[int]

    dv.solve(x_mat)
    print (read_mat.text_to_float('/home/peter/MISC/mat_file'))

    nrows = 4
    ncols = 4
    x_mat = util.generate_random_symmetric_matrix(nrows, ncols)
    print ('x_mat \n', x_mat )
    x_mat_scaled = util.make_diagonally_dominant(x_mat)
    print ('x_mat_scaled \n', x_mat_scaled)

main()


