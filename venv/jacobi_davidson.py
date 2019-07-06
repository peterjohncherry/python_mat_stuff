import numpy as np

def jacobi_davidson(A, num_eig, thresh, max_it ):

    ne = np.size(A,0)
    v_init = np.zeros(ne)
    v_init[0] = 1.0
    print ("v_init = " , v_init )

    #set initial number of guess vectors
    if 2*num_eig < ne:
        ngvecs = 2*num_eig
    else:
        ngvecs = num_eig

    v_vecs = np.eye(ne, ngvecs) # set of ngv unit vectors as guess
    w_vecs = np.zeros((ne, ngvecs))  # set of ngv unit vectors as guess

    for ii in range(ngvecs):
        w_vecs[:,ii] = A.dot(v_vecs[:,ii])

    for ii in range(ngvecs):
        print ("v_vecs[:,",ii,"]", v_vecs[:,ii] )
        print("w_vecs[:,", ii, "]", w_vecs[:, ii])


    V = np.zeros((ne,ne))# array of zeros to hold guess
    I = np.eye(ne) #identity matrix of same dimension as A
