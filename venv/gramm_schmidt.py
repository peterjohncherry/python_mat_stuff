import numpy as np

#THESE ALL ASSUME VECTORS ARE IN COLUMNS!!!!

# proj_{u}(v) = [(v|u)/(u|u)]*u
def proj_op(u,v):
    return ((v.dot(u))/ (u.dot(u)))*u

###########################################################################
#Unmodified Gramm-Schmidt algorithms
###########################################################################
#orthogonalizes u with respect to V
def gs_one_vec(V,u):
    new_u = u
    print("u = " , u )
    for ii in range(np.size(V,1)):
        new_u -= proj_op(u,V[:,ii])
        print("new_u = ", new_u)
    return new_u

def gs_one_vec_bad_v(V,u):
    new_u = u
    X = V
    gs_full(X)
    gs_one_vec(X, u)
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
        U[:,kk] = tmp
    return U

###########################################################################
#Modified Gramm-Schmidt algorithms
###########################################################################