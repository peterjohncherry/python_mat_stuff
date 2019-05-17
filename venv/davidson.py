import numpy as np
import cmath
import mat_utils as util

#crappy guess vector obtained by slicing row and normalizing
def get_guess_vec(mat_in, ii):
    nrows = mat_in.shape[0]
    return util.normalize(mat_in[ii,0:]).reshape(nrows, 1)

# Davidson algorithm
# AA is the matrix whose eigenvectors we seek
# uu is the list of initial guess vectors
# ww is the list of vectors defined by AA.uu[j]
def solve(AA):
    uu = []
    ww = []
    BB = np.array([[0.0]])
    print ("AA.shape = ", AA.shape[0] )
    for jj in range(AA.shape[0]):
        print (jj)
        uu.append(get_guess_vec(AA, jj))
        ww.append(util.multiply(AA, uu[jj]))
        print "jj = " +str(jj) + "\n"
        for kk in range(0, jj):
            BB[kk][jj] = util.dot(uu[kk], ww[jj])
            BB[jj][kk] = util.dot(uu[jj], ww[kk])
        BB_tmp = np.array(np.r_['-1', BB, np.zeros((BB.shape[1], 1))])
        BB = np.array(np.r_['0', BB_tmp, np.zeros((1, BB_tmp.shape[1]))])


