import numpy as np
import cmath
import davidson
import mat_utils as util
import plot_tools as pt
import matplotlib.pyplot as plt


def main():
    # define dimensions of the array
    nrows = 5
    ncols = 5
    x_mat = util.generate_random_matrix(nrows, ncols)  # type: List[int]
    print ("x")
    print (x_mat)

    y_vec = util.generate_random_vector(ncols)
    print "\n y \n"
    print (y_vec)

    z_vec = util.multiply(x_mat, y_vec)
    print "\n z \n"
    print (z_vec)

    q_scalar = util.multiply(y_vec,z_vec)
    print "\n q "
    print (q_scalar)

    mat_plots = pt.Matrix_Images('x_mat')
    mat_plots.add_plot(x_mat)
    x_mat_symm = util.symmetrize_mat(x_mat)
    mat_plots.add_plot(x_mat_symm)
    mat_plots.save_imgs()



main()

