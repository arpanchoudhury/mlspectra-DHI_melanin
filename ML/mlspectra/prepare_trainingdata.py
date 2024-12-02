import numpy as np
import kernels

def prepare_trainingdata(kernel,N_train,load_K,file_kernel,indices,lamd,X,Prop,opt_sigma):

    if load_K:
        K = np.load(file_kernel)
    else:
        K = np.zeros([N_train,N_train], dtype=float)
        print('Calculating kernel elements'+'\n')
        if kernel == 'laplacian':
            for itrain in range(N_train):
                K[itrain,itrain] = 1.0 + lamd

                for jtrain in range(itrain+1,N_train):
                    Xt = X[indices[itrain]]
                    Xq = X[indices[jtrain]]
                    Yt = np.zeros([1,len(Xt)],dtype=float)
                    Yq = np.zeros([1,len(Xq)],dtype=float)
                    Yt[0] = Xt
                    Yq[0] = Xq
                    tmp = kernels.laplacian_kernel(Yq,Yt, sigma=opt_sigma)
                    K[itrain,jtrain] = tmp[0,0]
                    K[jtrain,itrain] = K[itrain,jtrain]
            np.save(file_kernel, K)
        elif kernel == 'gaussian':
            for itrain in range(N_train):
                K[itrain,itrain] = 1.0 + lamd

                for jtrain in range(itrain+1,N_train):
                    Xt = X[indices[itrain]]
                    Xq = X[indices[jtrain]]
                    Yt = np.zeros([1,len(Xt)],dtype=float)
                    Yq = np.zeros([1,len(Xq)],dtype=float)
                    Yt[0] = Xt
                    Yq[0] = Xq
                    tmp = kernels.gaussian_kernel(Yq,Yt, sigma=opt_sigma)
                    K[itrain,jtrain] = tmp[0,0]
                    K[jtrain,itrain] = K[itrain,jtrain]
            np.save(file_kernel, K)

        elif kernel == 'linear':
            for itrain in range(N_train):
                K[itrain,itrain] = 1.0 + lamd
                if np.mod(itrain,10) == 0:
                    print(itrain, 'rows calculated', N_train-itrain, 'remaining')
                for jtrain in range(itrain+1,N_train):
                    Xt = X[indices[itrain]]
                    Xq = X[indices[jtrain]]
                    Yt = np.zeros([1,len(Xt)],dtype=float)
                    Yq = np.zeros([1,len(Xq)],dtype=float)
                    Yt[0] = Xt
                    Yq[0] = Xq
                    tmp = kernels.linear_kernel(Yq,Yt, sigma=opt_sigma)
                    K[itrain,jtrain] = tmp[0,0]
                    K[jtrain,itrain] = K[itrain,jtrain]
            np.save(file_kernel, K)



    N_prop = len(Prop[0])
    P = np.zeros([N_train,N_prop],dtype=float)

    for iprop in range(N_prop):
        for itrain in range(N_train):
            P[itrain,iprop] = Prop[indices[itrain],iprop]


    return K, P

