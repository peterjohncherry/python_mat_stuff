import numpy as np

def proj_op(v,u):
    return (v.dot(u)/u.dot(u))*u


#orthogonalizes u with respect to V
def gs_one_vec(V, u):
    new_u = u
    for ii in range(len(u)):
        new_u -= proj_op(V[:,ii],u)
    return new_u


