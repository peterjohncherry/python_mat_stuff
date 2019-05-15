import numpy as np
import cmath
import davidson
import mat_utils as util


def main():
    mat = util.generate_random_symmetric_matrix(5,5)
    new_mat = util.symmetrize_mat(mat)
    eigvals, eigvecs =  util.diag( new_mat )

    print (mat)
    print "\n"
    x = util.generate_random_matrix(5,1)  # type: List[int]
    print (x)
    print "\n"
    y = np.matmul(mat, x)
    print (y)

main()

