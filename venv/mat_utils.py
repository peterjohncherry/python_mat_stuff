import numpy as np
import cmath


# generates a real random matrix with nrows and ncols
def generate_random_matrix(nrows, ncols):
    # type: (int, int) -> matrix
    return np.random.rand(nrows, ncols)


# symmetrizes matrix so it diagonalizes nice like
def symmetrize_mat(mat):
    # type : (numpy.matrix) -> numpy.matrix
    new_mat = mat / 2 + mat.transpose() / 2
    return new_mat


# generates a real random matrix with nrows and ncols
def generate_random_symmetric_matrix(nrows, ncols):
    # type: (int, int) -> matrix
    return symmetrize_mat(np.random.rand(nrows, ncols))


# diagonalizes a matrix directly
def diag(mat):
    return np.linalg.eig(mat)

def multiply(mat_in, vec_in):
    return np.mat_mult(mat_in, vec_in)

