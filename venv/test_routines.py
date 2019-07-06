import davidson as dv
import jacobi_davidson as jd
import mat_utils as util
import numpy as np

def test_davidson():
    tol = 1e-8  # Convergence tolerance
    maxit = 30

    ne = 1200
    sparsity = 0.000001
    A = np.zeros((ne, ne))
    for i in range(0, ne):
        A[i, i] = i + 1

    A = A + sparsity * np.random.randn(ne, ne)
    A = (A.T + A) / 2
    eig = 4  # number of eigenvalues to solve
    dv.davidson(A, tol, maxit, eig)

def test_jacobi_davidson():
    tol = 1e-8  # Convergence tolerance
    maxit = 30

    ne = 10
    sparsity = 0.000001
    A = np.zeros((ne, ne))
    for i in range(0, ne):
        A[i, i] = i + 1

    A = A + sparsity * np.random.randn(ne, ne)
    A = (A.T + A) / 2
    eig = 4  # number of eigenvalues to solve
    jd.jacobi_davidson(A, eig, tol, maxit)

def test_reader():
    print (read_mat.text_to_float('/home/peter/MISC/mat_file'))

def test_math_generation():
    nrows = 4
    ncols = 4
    x_mat = util.generate_random_matrix(nrows, ncols)  # type: List[int]
    dv.solve(x_mat)
    x_mat = util.generate_random_symmetric_matrix(nrows, ncols)
    print ('x_mat \n', x_mat )
    x_mat_scaled = util.make_diagonally_dominant(x_mat)
    print ('x_mat_scaled \n', x_mat_scaled)

def test_solver():
    x_mat_scaled = util.make_diagonally_dominant(util.generate_random_symmetric_matrix(nrows, ncols))
    print('x_mat_scaled \n', x_mat_scaled)
    dv.solve(x_mat_scaled)