import numpy as np
import cmath
import mat_utils as util

#crappy guess vector obtained by slicing row and normalizing
def get_guess_vec(mat_in, ii):
    nrows = mat_in.shape[0]
    return util.normalize(mat_in[ii,0:]).reshape(1, nrows)[0]

# Davidson algorithm
# AA is the matrix whose eigenvectors we seek
# uu is the list of initial guess vectors
# ww is the list of vectors defined by AA.uu[j]
def solve(AA):
    uu = np.array([get_guess_vec(AA,0)])
    print "uu", uu
    ww = np.array(util.multiply(uu, AA))

    print "ww", ww
    BB = np.array([ util.dot(uu, ww) ])                                 # b[j][j] = (u[j]|w[j])
    for jj in range(AA.shape[0]):
        if jj != 0 :
            uu = np.r_['0', uu, np.array([get_guess_vec(AA,jj)])]
            print ("uu = ", uu)
            print ""
            ww = util.multiply(AA, uu.transpose())                     #w[j] = A * u[j]
            print ("ww = ", ww)
            BB = util.multiply(ww, uu)
            print "BB = ", BB
            eigvals, eigvecs = util.diag(BB)
            print eigvecs, eigvals
            ss = eigvecs[np.argmax(eigvals),:].transpose()
            print "ss = ", ss