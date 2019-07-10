import numpy as np
import gramm_schmidt as gs
import mat_utils as util

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
    print ("v_vecs = ", v_vecs)
    w_vecs = np.zeros((ne, 1)) # w[i] = Av[i]

    for ii in range(1):
        w_vecs[:,ii] = A.dot(v_vecs[:,ii])

    H = build_H(v_vecs, w_vecs)
    thetas = np.ndarray((ngvecs))
    thetas[0] = H[0,0]
    ritz_vec = v_vecs[:,0]
    residual_vec = w_vecs[:,0] - thetas[0]*ritz_vec
    residual_norm = np.linalg.norm(residual_vec)

    approx_type = "diag(A)"
    iter = 1

    #2. Iterate until convergence
    while iter < max_it and residual_norm > thresh:
        print  ( "\n\n iter = ", iter )
        #3. Inner Loop
        M = build_M(A, ritz_vec, thetas, approx_type)

        t = np.matmul(np.linalg.inv(M), -residual_vec)
        t = gs.modified_gs_one_vec(v_vecs, t)
        t = np.reshape(t,(ne,1))

        v_vecs = np.append(v_vecs, t, axis=1)
        v_vecs.reshape((ne, iter+1))
        v_vecs = gs.modified_gs_full(v_vecs)
        util.check_column_orthogonality(v_vecs, thresh=0.000000001, name=" V")

        w_vecs = np.append(w_vecs, np.matmul(A,t), axis =1)
        w_vecs.reshape((ne, iter+1))

        H = np.matmul(v_vecs.transpose(), w_vecs)

        # get eig vals and eigvecs of new H and sort
        thetas, s_vecs = np.linalg.eig(H)
        idx = thetas.argsort()
        thetas = thetas[idx]
        s_vecs = s_vecs[:, idx]
        print("thetas = ", thetas)
        print ("s_vecs = ", s_vecs)

        ritz_vec = np.matmul(v_vecs, s_vecs[:,0])
        uhat_vec = np.matmul(A, ritz_vec)
        residual_vec = uhat_vec-thetas[0]*ritz_vec

        residual_norm = np.linalg.norm(residual_vec)

        if residual_norm > thresh :
            print ("||r||> thresh  :: ", residual_norm, ">", thresh, " -----> Not yet converged")
           # v_vecs[:,0] = ritz_vec
           # w_vecs[:,0] = uhat_vec
           # v_vecs = gs.modified_gs_full(v_vecs)
           # w_vecs = gs.modified_gs_full(w_vecs)
        else:
            print("||r|| <= thresh  :: ", residual_norm, "<=", thresh, " -----> Converged!")
            print ("thetas = ", thetas)
            break
        iter += 1





