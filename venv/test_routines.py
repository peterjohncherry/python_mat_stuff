import davidson as dv
import mat_utils as util
import numpy as np

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