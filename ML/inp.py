import numpy as np
from sklearn.model_selection import train_test_split
import mlspectra
import random



#---------------- Data loading ------------------------------

X = np.loadtxt("fingerprint_64.dat")
y = np.loadtxt("spectra_12bins_waveuniform.dat")
y = y[:-2,:] # last two rows contain bin boundary and binwidths


n_times_list = []
for N_times in range(20):
    X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, test_size=210)
    indices = list(i_mol for i_mol in range(X_train.shape[0]))
    indices_q = list(i_mol for i_mol in range(X_test.shape[0]))


    #--------------- Kernel specific inputs ---------------------------
    kernel = 'laplacian'
    load_K = False
    file_kernel = 'kernel.npy'
    opt_sigma = 10.0
    lamd = 1e5

    #------------------- Training & prediction -------------------------
    obs_list  = [] 
    for N_train in [1,10,100,500,1000,5000,10000]:
        K, P = mlspectra.prepare_trainingdata(kernel,N_train,load_K,file_kernel,indices,lamd,X_train,y_train,opt_sigma) 
        N_prop = y_train.shape[1]

        print('solving matrix equation...')
        alpha = mlspectra.linalg_solve(N_prop,K,P,solver='inv')
    
        print('predicting...')
        out_of_sample_obs = []
	out_of_sample_mae = []

        y_pred_arr = []
        for iquery in range(X_test.shape[0]):
            y_pred, _ = mlspectra.predict(kernel,X_train,X_test,alpha,indices,indices_q,iquery,opt_sigma)
            y_act = y_test[indices_q[iquery],:]

	    print(y_pred.shape[0])
            for j in range(y_pred.shape[0]):
                if y_pred[j] < 0:
                    y_pred[j] = 0.0
                

	    mae = np.abs(y_pred - y_act)
	    y_pred_arr.append(y_pred)

            y_act_sum = np.linalg.norm(y_act)
            y_pred_sum = np.linalg.norm(y_pred)
	    print('act',y_act)
            print('pred',y_pred)

            print('act sum',y_act_sum)
            print('pred sum',y_pred_sum)

            y_act_norm = y_act/y_act_sum
            y_pred_norm = y_pred/y_pred_sum

            phi = np.linalg.norm(np.abs(y_pred_norm - y_act_norm))
            obs = 1 - phi
            out_of_sample_obs.append(obs)
            out_of_sample_mae.append(mae)

        #y_pred_arr = np.array(y_pred_arr)
        #np.savetxt("y_pred_12bins_waveuniform_X-fing40_linearkernel.dat", y_pred_arr)


        mean_obs = np.mean(out_of_sample_obs)
	mean_mae = np.mean(out_of_sample_mae, axis=0)
	
	with open(kernel+'_fingerprint_64_12bins_waveuniform_ntimes_withcrossMAEHyperparams_MAE.dat', "a+") as res:
            res.write(str(mean_mae))
        obs_list.append(mean_obs)
    n_times_list.append(obs_list)

n_times_list = np.array(n_times_list)
np.savetxt(kernel + '_OBS_fingerprint_64_12bins_waveuniform_ntimes_withcrossMAEHyperparams_norm.dat',n_times_list)
