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
        uu.append(get_guess_vec(AA, jj))
        ww.append(util.multiply(AA, uu[jj]))            #w[j] = A * u[j]
        for kk in range(0, jj):
            BB[kk][jj] = util.dot(uu[kk], ww[jj])                           #b[k][j] = (u[k]|w[j])
            BB[jj][kk] = util.dot(uu[jj], ww[kk])                           #b[j][k] = (u[j]|w[k])
        BB[jj][jj] = util.dot(uu[jj], ww[jj])                               #b[j][j] = (u[j]|w[j])
        BB_tmp = np.array(np.r_['-1', BB, np.zeros((BB.shape[1], 1))])
        BB = np.array(np.r_['0', BB_tmp, np.zeros((1, BB_tmp.shape[1]))])
        eigvals, eigvecs = util.diag(BB)
        print ( "eigvals = ", eigvals)
        print ("eigvecs = ", eigvecs )
        SS = util.normalize(eigvecs[:, np.argmax(eigvals)])
        #YY = util.multiply(uu,SS)        # Ritz vector
        #AAYY = util.multiply(AA, YY)     #
        #residual = AAYY- eigvals[np.argmax(eigvals)]*YY
        #TT = np.diagonal(a).copy()
        #for ii in range(jj):
        #    TT[ii] = eigvals[ii]
        #    TT[ii] = 1/TT[ii]


