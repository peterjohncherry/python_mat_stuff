import numpy as np
import cmath
import davidson
import mat_utils as util



def main():
    print ( "hello")
    mat = util.generate_random_symmetric_matrix(3,3)
    print mat
    new_mat = util.symmetrize_mat(mat)
    print new_mat
    eigvals, eigvecs =  util.diag( new_mat )
    print ("eigvecs\n")
    print (eigvecs)
    print ("eigvals\n")
    print (eigvals)

main()

