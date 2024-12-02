import numpy as np
import kernels


def predict(kernel,X_train, X_query,alpha,indices_t,indices_q,iquery,opt_sigma):
    N_train = alpha.shape[0]
    N_prop = alpha.shape[1]

    K = np.zeros(N_train)

    if kernel == 'laplacian':
        for itrain in range(N_train):
            Xt = X_train[indices_t[itrain]]
            Xq = X_query[indices_q[iquery]]
            Yt = np.zeros([1,len(Xt)],dtype=float)
            Yq = np.zeros([1,len(Xq)],dtype=float)
            Yt[0] = Xt
            Yq[0] = Xq
            tmp = kernels.laplacian_kernel(Yq,Yt, sigma=opt_sigma)
            K[itrain] = tmp[0,0]
    elif kernel == 'gaussian':
        for itrain in range(N_train):
            Xt = X_train[indices_t[itrain]]
            Xq = X_query[indices_q[iquery]]
            Yt = np.zeros([1,len(Xt)],dtype=float)
            Yq = np.zeros([1,len(Xq)],dtype=float)
            Yt[0] = Xt
            Yq[0] = Xq
            tmp = kernels.gaussian_kernel(Yq,Yt, sigma=opt_sigma)
            K[itrain] = tmp[0,0]
    elif kernel == 'linear':
        for itrain in range(N_train):
            Xt = X_train[indices_t[itrain]]
            Xq = X_query[indices_q[iquery]]
            Yt = np.zeros([1,len(Xt)],dtype=float)
            Yq = np.zeros([1,len(Xq)],dtype=float)
            Yt[0] = Xt
            Yq[0] = Xq
            tmp = kernels.linear_kernel(Yq,Yt, sigma=opt_sigma)
            K[itrain] = tmp[0,0]


    P_pred = np.zeros(N_prop)
    for i_prop in range(N_prop):
        P_pred[i_prop] = np.dot(K,alpha[:,i_prop])

    return P_pred, Xq
