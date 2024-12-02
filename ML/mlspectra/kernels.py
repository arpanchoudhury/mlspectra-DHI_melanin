import numpy as np
import fortran_kernels

def laplacian_kernel(A, B, sigma):

    na = A.shape[0]
    nb = B.shape[0]

    K = np.empty((na, nb), order='F')

    fortran_kernels.flaplacian_kernel(A.T, na, B.T, nb, K, sigma)

    return K

def gaussian_kernel(A, B, sigma):

    na = A.shape[0]
    nb = B.shape[0]

    K = np.empty((na, nb), order='F')

    fortran_kernels.fgaussian_kernel(A.T, na, B.T, nb, K, sigma)

    return K

def linear_kernel(A, B, sigma):

    na = A.shape[0]
    nb = B.shape[0]

    K = np.empty((na, nb), order='F')

    fortran_kernels.flinear_kernel(A.T, na, B.T, nb, K)

    return K
