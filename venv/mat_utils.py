import numpy as np
import cmath

#generate random array
def generate_random_vector( vec_length):
    return np.random.random_sample((vec_length,))

# generates a real random matrix with nrows and ncols
def generate_random_matrix(nrows, ncols):
    # type: (int, int) -> matrix
    return np.random.rand(nrows, ncols)

# symmetrizes matrix
def symmetrize_mat(mat):
    # type : (numpy.matrix) -> numpy.matrix
    new_mat = mat / 2 + mat.transpose() / 2
    return new_mat

# generates a real random matrix with nrows and ncols
def generate_random_symmetric_matrix(nrows, ncols):
    # type: (int, int) -> matrix
    return symmetrize_mat(np.random.rand(nrows, ncols))

#Scales down off diagonal elements to make matrix amenable to solution by Davidson
def make_diagonally_dominant(mat_in) :
    mat_out = mat_in
    for ii in range(np.size(mat_in,0)):
        for jj in range(np.size(mat_in, 1)):
            if ii != jj:
                mat_out[ii,jj] = mat_in[ii,jj]/(1000**abs(ii-jj))
    return mat_out

# diagonalizes a matrix directly
def diag(mat):
    return np.linalg.eig(mat)

# matrix multiplication, will work either with matrices or matrix vectors
def multiply(mat_in, vec_in):
    return np.matmul(mat_in, vec_in)

# dot product
def dot(vec_in1, vec_in2):
    return np.vdot(vec_in1, vec_in2)

#get the Frobenius norm
def norm(array_in):
    return np.linalg.norm(array_in)

#get sub matrix, note zero indexing
def get_sub_matrix( mat_in, top_row, bottom_row, left_column, right_column):
    return mat_in[np.ix_(range(top_row, bottom_row+1), range(left_column, right_column+1))]

#get sub matrix, takes ranges as inputs, note zero indexing
def get_sub_matrix( mat_in, row_range, col_range):
    return mat_in[np.ix_(row_range, col_range)]

#normalizes vector, but will return original vector if norm is 0
def normalize( vec_in ):
    my_norm = norm(vec_in)
    if my_norm == 0:
        return vec_in
    return vec_in/(my_norm)

#checks whether the columns of A are orthogonal to one another
def check_column_orthogonality(A, thresh = 0.000000001, name="???"):
    ncols = np.size(A,1)
    for jj in range(ncols):
        for ii in range(ncols-1):
            if ii != jj :
                self_dot =  np.abs(np.dot(A[:,jj],A[:,ii]))
                if self_dot > thresh :
                    print("A[:,",jj,"] dot A[:,",ii,"] = ", self_dot)
                    print ("ORTHOGONALIZATION OF", name, "FAILED!!")
                    return

#checks whether the columns of A are orthogonal to one another
def check_column_orthogonality_u(A, u, thresh = 0.000000001, name="???"):
    ncols = np.size(A,1)
    for jj in range(ncols):
        self_dot = np.abs(np.dot(A[:,jj],u))
        if self_dot > thresh :
            print("u dot A[:,",jj,"] = ", self_dot)
            print ("ORTHOGONALIZATION OF", name, "FAILED!!")
            return
