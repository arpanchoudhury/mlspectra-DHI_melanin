import numpy as np
from scipy import linalg
from scipy.linalg import cho_factor, cho_solve

def linalg_solve(N_prop,K,P,solver='cholesky'):
    print(K.shape)
    N_train=K.shape[0]

    alpha=np.zeros([N_train,N_prop])

    if solver == 'cholesky':
        Klow, low = cho_factor(K)
        print('Klow, low done')
        print('solving cho solve..')
        for i_prop in range(N_prop):
            alpha[:,i_prop] = cho_solve((Klow, low), P[:,i_prop])
    else:
        alpha = linalg.solve(K, P)

    return alpha
