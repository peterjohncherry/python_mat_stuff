import numpy as np

def jacobi_davidson(A, num_eig, thresh, max_it ):

    ne = np.size(A,0)
    v_init = np.zeros(ne)
    print ("v_init = " , v_init )
