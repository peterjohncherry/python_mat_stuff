import numpy as np
import mat_utils as util

#THESE ALL ASSUME VECTORS ARE IN COLUMNS!!!!

# proj_{u}(v) = [(v|u)/(u|u)]*u
def proj_op(u,v):
    return ((v.dot(u))/(u.dot(u)))*u

###########################################################################
#Unmodified Gramm-Schmidt algorithms
###########################################################################
#orthogonalizes u with respect to V
def gs_one_vec(V,u):
    new_u = u
    ncols = np.size(V,1)
    for ii in range(ncols):
        new_u = new_u - proj_op(V[:, ii], new_u)
    return new_u

def gs_one_vec_bad_v(V,u):
    X = V
    gs_full(X)
    util.check_column_orthogonality(X, name = "bob")
    new_u = gs_one_vec(X, u)
    return new_u

#orthogonalizes u with respect to V
def gs_full(V):
    U = V
    ncols = np.size(V,1)
    for kk in range(1, ncols):
        tmp = V[:, kk]
        jj = 0
        while jj < kk:
            tmp = tmp - proj_op(U[:, jj], V[:, kk])
            jj += 1
        U[:, kk] = tmp
    return U

###########################################################################
#Modified Gramm-Schmidt algorithms
###########################################################################

def modified_gs_full(X):
    V = X
    Q = np.zeros_like(V)
    ncols = np.size(V, 1)
    ii = 0
    while ii < ncols :
        rii = np.sqrt(np.dot(V[:,ii], V[:,ii]))
        Q[:,ii] = V[:,ii]/rii
        jj = ii
        while jj < ncols:
            rij = np.dot(Q[:,ii], V[:,jj])
            V[:,jj] = V[:,jj]-rij*Q[:,ii]
            jj +=1
        ii+=1
    return Q

def modified_gs_one_vec(X,u):
    V = X
    ncols = np.size(V, 1)
    rii = np.sqrt(np.dot(u, u))
    q = u/rii
    jj = 0
    while jj < ncols:
        rij = np.dot(q, V[:,jj])
        q = q - rij*V[:,jj]
        jj +=1
    return q