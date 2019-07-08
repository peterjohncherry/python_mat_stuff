import numpy as np
import gramm_schmidt as gs

def build_H(v_vecs, w_vecs):
    ne = np.size(v_vecs, 1)
    H = np.ndarray((ne, ne))
    for ii in range(np.size(v_vecs, 1)):
        for jj in range(np.size(w_vecs, 1)):
            H[ii, jj] = v_vecs[:, ii].dot(w_vecs[:, jj])
    return H

#get diagonal matrix, DA, of size dim
def get_DA(A, dim):
    DA = np.ndarray((dim,dim))
    for ii in range(dim):
        DA[ii,ii] = A[ii,ii]
    return DA

# Approximation for (I-uu*)(A-theta*I)(I-uu*)
# Well, currently exact
def build_M( A, u_vecs, thetas, approx_type):
    ne = np.size(A, 1)
    if approx_type == "exact" :
        uu = np.matmul(u_vecs, np.transpose(u_vecs))
        I = np.eye(ne,ne)
        projector = I-uu

        tmp = np.ndarray((ne, ne))
        for ii in range(ne):
            for jj in range(ne):
                tmp[ii,jj] = A[ii, ii] - thetas[ii]

        return np.matmul(projector, np.matmul(tmp, projector))

    elif approx_type == "diag(A)":
        return get_DA(A, ne)

# The variable names and numbers in comments correspond to those in Algorithm 1 in
# Sleijpen & Van der Vorst, SIAM Review, Vol. 42, No. 2, 2000, 267-293
def jacobi_davidson(A, num_eig, thresh, max_it ):
    # 1. Start
    ne = np.size(A,0)

    #set initial number of guess vectors
    if 2*num_eig < ne:
        ngvecs = 2*num_eig
    else:
        ngvecs = num_eig

    v_vecs = np.eye(ne, 1) # set of ngv unit vectors as guess
    w_vecs = np.zeros((ne, 1)) # w[i] = Av[i]

    for ii in range(1):
        w_vecs[:,ii] = A.dot(v_vecs[:,ii])

    u_vecs = v_vecs

    H = build_H(v_vecs, w_vecs)
    thetas = np.ndarray((ngvecs))
    thetas[0] = H[0,0]

    r = w_vecs[:,0] - thetas[0]*u_vecs[:,0]
    diff = 2* thresh

    approx_type = "diag(A)"
    diff = 100000
    iter = 1
    max_it = 5
    #2. Iterate until convergence
    while iter < max_it  and diff > thresh  :
        #diff = thresh
        #3. Inner Loop
        M = build_M(A, u_vecs, thetas, approx_type)
        Minv = np.linalg.inv(M)
        t = np.matmul(np.linalg.inv(M), r)
        t = gs.gs_one_vec(v_vecs, t)
        t = np.reshape(t,(ne,1))
        if iter == 1:
            v_vecs = t
            w_vecs = np.matmul(A,t)
        else:
            v_vecs = np.append(v_vecs, t,axis=1)
            v_vecs.reshape((ne,iter))
            w_vecs = np.append(w_vecs, np.matmul(A,t), axis =1 )
            v_vecs.reshape((ne, iter))
        iter += 1
        print ("w_vecs = \n", w_vecs)

        diff += 1

       # w[:,kk]
