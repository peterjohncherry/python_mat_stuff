import numpy as np
import cmath
import mat_utils as util

#crappy guess vector obtained by slicing row and normalizing
def get_guess_vec(mat_in, ii):
    nrows = mat_in.shape[0]
    print ("nrows in guess vec = " + str(nrows) + "\n")
    return util.normalize(mat_in[ii,0:]).reshape(nrows, 1)

# Davidson algorithm
# AA is the matrix whose eigenvectors we seek
# uu is the list of initial guess vectors
# ww is the list of vectors defined by AA.uu[j]
def solve(AA):
    uu = []
    ww = []
    for ii in range(AA.shape[0]):
        uu.append(get_guess_vec(AA,ii))
        ww.append(util.multiply(AA, uu))
        print " uu["+str(ii) + "] \n", uu[ii], "\n"
        print " ww["+str(ii) + "] \n", ww[ii], "\n"



