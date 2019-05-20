import numpy as np
import cmath
import davidson as dv
import mat_utils as util
import plot_tools as pt
import matplotlib.pyplot as plt


def main():
    # define dimensions of the array
    nrows = 6
    ncols = 6
    x_mat = util.generate_random_matrix(nrows, ncols)  # type: List[int]

    dv.solve(x_mat)

main()

