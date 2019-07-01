import numpy as np

def my_goings():
    #ne = 1200  # Dimension of matrix
    #tol = 1e-8  # Convergence tolerance
    maxit = 30

    ne = 4
    sparsity  = 0.000001
    A = np.zeros((ne,ne))
    for i in range (0,ne):
        A[i,i] = i+1
    A = A+sparsity*np.random.randn(ne,ne)
    A = (A.T +A)/2

    ngv =8 # number of guess vectors
    eig = 4 # number of eigenvalues to solve
    gvecs = np.eye(ne,ngvecs) # set of ngv unit vectors as guess
    V = np.zeros((ne,ne))# array of zeros to hold guess
    I = np.eye(ne) #identity matrix of same dimension as A


    x = range(ngvecs, maxit)
    for mm in range ( ngvecs, maxit, ngvecs):
        if mm == ngvecs : #check if first iteration
            for jj in range (0,ngvecs): # build up initial guess vectors
                V[:, j] = gvecs[:,jj]/np.linalg.norm(t[:,jj])
            theta_old = 1 # arbitrary value for initial eigenvalue for comparison

        elif mm > ngvecs : # if not first iteration, set theta_old to eigvals from last iteration
            theta_old = theta[:eig]

        # Compute QR factorization V is matrix with orthonormal columns
        #                          R is upper triangular matrix
        V, R = np.linalg.qr(V)

        # project matrix A onto subspace defined by new guess vectors V
        VT = V[:,:(m+1)].T
        AV = np.(A,V[:,:(m+1)])
        T = np.dot(VT, AV)

        THETA, S =`

