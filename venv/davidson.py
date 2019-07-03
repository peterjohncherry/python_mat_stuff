import numpy as np

def davidson():

    tol = 1e-8  # Convergence tolerance
    maxit = 30

    ne = 1200
    sparsity = 0.000001
    A = np.zeros((ne,ne))
    for i in range (0,ne):
        A[i,i] = i+1

    A = A+sparsity*np.random.randn(ne,ne)
    A = (A.T +A)/2

    ngvecs =8 # number of guess vectors
    eig = 4 # number of eigenvalues to solve
    gvecs = np.eye(ne,ngvecs) # set of ngv unit vectors as guess
    V = np.zeros((ne,ne))# array of zeros to hold guess
    I = np.eye(ne) #identity matrix of same dimension as A


    x = range(ngvecs, maxit)
    for mm in range ( ngvecs, maxit, ngvecs):
        if mm == ngvecs :                # check if first iteration
            for jj in range (0,ngvecs):  # build up initial guess vectors
                V[ :, jj] = gvecs[:,jj]/np.linalg.norm(gvecs[:,jj])
            theta_old = 1                # arbitrary value for initial eigenvalue for comparison

        elif mm > ngvecs : # if not first iteration, set theta_old to eigvals from last iteration
            theta_old = theta[:eig]

        # Compute QR factorization V is matrix with orthonormal columns
        #                          R is upper triangular matrix
        V, R = np.linalg.qr(V)

        # project matrix A onto subspace defined by new guess vectors V
        VT = V[ :, :(mm+1)].T
        AV = np.dot( A , V[ :, :(mm+1)])
        VTAV = np.dot(VT, AV)

        #Get eigenvectors and eigenvalues and sort them
        THETA, S = np.linalg.eig(VTAV)
        idx = THETA.argsort()
        theta = THETA[idx]
        s = S[:, idx]

        #  calculate residual and new vector, Extend search space, check convergence
        for jj in range (0,ngvecs):
            uj = np.dot( V[ :, :(mm+1) ], s[ :, jj] ) # V.sj
            AmtI =A-theta[jj]*I                       # (A-I.theta)
            rj = np.dot(AmtI, uj)                     #  rj = (A.uj -theta.uj)
            tj = rj/(theta[jj]- A[jj,jj])             # get t from (DA.uj -theta.uj)tj = rj ; DA is diagonal of A
            V[:, (mm+jj+1)] = tj                      # Add this tj to search space
        norm = np.linalg.norm(theta[:eig]-theta_old)  # Calculate change from eigvals on last iteration
        if norm < tol:                                # Check convergence
            break

    print ("davidson results = ", theta[:eig])

    E, Vec = np.linalg.eig(A)
    E = np.sort(E)
    print ("numpy results = ", E[:eig])

   # mat_sections()