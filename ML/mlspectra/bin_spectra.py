import pandas as pd
import numpy as np
import mlspectra

def bin_spectra_uniform(spec_path, read_P, file_spec, file_geom, e_min, e_max, N_bin, N_fing, spec_files, dihedral_file):
    #spec_files = mlspectra.read_files(spec_path)
    fo = open(spec_files, 'r')
    filenames = fo.readlines()

    dihedrals = np.loadtxt(dihedral_file)
    dihedrals = np.abs(dihedrals)

    #e_max = 6.1992
    #e_min = 1.5498
    #N_bin=30

    lambda_min=[]
    lambda_max=[]
    dlambda = (e_max - e_min)/N_bin
    print(dlambda)

    for i_bin in range(N_bin):
        lambda_min.append( e_min + (i_bin)*dlambda )
        lambda_max.append( e_min + (i_bin+1)*dlambda )

    for i in range(len(lambda_min)):
        lambda_min[i] = 4.1357*2.9979*1e2/lambda_min[i] 
        lambda_max[i] = 4.1357*2.9979*1e2/lambda_max[i]

    lambda_min = lambda_min[-1::-1]
    lambda_max = lambda_max[-1::-1]

    for i in range(len(lambda_min)):
        print(lambda_max[i], lambda_min[i])

    width = []
    for i in range(len(lambda_min)):
        width.append(lambda_min[i]-lambda_max[i])

    #plt.bar(lambda_max, 1.0, width=width, bottom=None, align='edge', edgecolor='black', color='white')
    #plt.show()

    N_file = len(filenames)
    i_file = 0

    Int_lam=np.zeros([N_file+2,N_bin])
    geom_mat = np.zeros([N_file,N_fing])
    if N_fing == 40:
        discrete = 2
    elif N_fing == 46:
        discrete = 4
    elif N_fing == 64:
        discrete = 10
    else:
        discrete = 0


    for i in filenames:
        spec_csv = i.strip()

        if np.mod(i_file, 100) == 0:
            print(i_file, ' out of ', N_file, ' done')
        spec_data=pd.read_csv(spec_path+'/'+spec_csv, names=['wavelength_nm', 'osc_strength'])

        wavelength=np.array(spec_data['wavelength_nm'])
        f=np.array(spec_data['osc_strength'])
        for i_bin in range(N_bin):
            sum_f=0.0
            for i_state in range(f.shape[0]):
                if wavelength[i_state] > lambda_max[i_bin] and wavelength[i_state] <= lambda_min[i_bin]:
                    sum_f=sum_f+f[i_state]
            Int_lam[i_file,i_bin]=sum_f
	print(spec_csv)
	pos_arr, ox_arr, pattern_arr, dihedral_arr = mlspectra.calfingerprint(spec_csv, dihedrals[i_file], discrete)
        #fprint = np.array([1,1,1,1])
        if N_fing == 32: 
            fprint = np.concatenate((pos_arr,ox_arr))
        elif N_fing == 34:
            fprint = np.concatenate((pos_arr,ox_arr,pattern_arr))
        elif N_fing == 38:
            fprint = np.concatenate((pos_arr,ox_arr,dihedral_arr))
        elif N_fing >= 40:
            fprint = np.concatenate((pos_arr,ox_arr,pattern_arr,dihedral_arr))

        geom_mat[i_file] = fprint
        i_file=i_file+1

    Int_lam[N_file,:] = lambda_min
    Int_lam[N_file+1,:] = dlambda
    np.savetxt(file_spec, Int_lam)    
    np.savetxt(file_geom, geom_mat)
    return Int_lam, lambda_min, dlambda



def bin_spectra_waveuniform(spec_path, read_P, file_spec, file_geom, e_min, e_max, N_bin, N_fing, spec_files, dihedral_file):
    #spec_files = mlspectra.read_files(spec_path)
    fo = open(spec_files, 'r')
    filenames = fo.readlines()

    dihedrals = np.loadtxt(dihedral_file)
    dihedrals = np.abs(dihedrals)

    #e_max = 6.1992
    #e_min = 1.5498
    #N_bin=30

    lambda_min=[]
    lambda_max=[]
    dlambda = (e_max - e_min)/N_bin
    print(dlambda)

    for i_bin in range(N_bin):
        lambda_min.append( e_min + (i_bin)*dlambda )
        lambda_max.append( e_min + (i_bin+1)*dlambda )

    #for i in range(len(lambda_min)):
    #    lambda_min[i] = 4.1357*2.9979*1e2/lambda_min[i]
    #    lambda_max[i] = 4.1357*2.9979*1e2/lambda_max[i]

    #lambda_min = lambda_min[-1::-1]
    #lambda_max = lambda_max[-1::-1]

    for i in range(len(lambda_min)):
        print(lambda_min[i], lambda_max[i])

    width = []
    for i in range(len(lambda_min)):
        width.append(lambda_min[i]-lambda_max[i])

    #plt.bar(lambda_max, 1.0, width=width, bottom=None, align='edge', edgecolor='black', color='white')
    #plt.show()

    N_file = len(filenames)
    i_file = 0

    Int_lam=np.zeros([N_file+2,N_bin])
    geom_mat = np.zeros([N_file,N_fing])
    if N_fing == 40:
        discrete = 2
    elif N_fing == 46:
        discrete = 4
    elif N_fing == 64:
        discrete = 10
    else:
        discrete = 0


    for i in filenames:
        spec_csv = i.strip()

        if np.mod(i_file, 100) == 0:
            print(i_file, ' out of ', N_file, ' done')
        spec_data=pd.read_csv(spec_path+'/'+spec_csv, names=['wavelength_nm', 'osc_strength'])

        wavelength=np.array(spec_data['wavelength_nm'])
        f=np.array(spec_data['osc_strength'])
        for i_bin in range(N_bin):
            sum_f=0.0
            for i_state in range(f.shape[0]):
                if wavelength[i_state] > lambda_min[i_bin] and wavelength[i_state] <= lambda_max[i_bin]:
                    sum_f=sum_f+f[i_state]
            Int_lam[i_file,i_bin]=sum_f
        print(spec_csv)
        pos_arr, ox_arr, pattern_arr, dihedral_arr = mlspectra.calfingerprint(spec_csv, dihedrals[i_file], discrete)
        #fprint = np.array([1,1,1,1])
        if N_fing == 32:
            fprint = np.concatenate((pos_arr,ox_arr))
        elif N_fing == 34:
            fprint = np.concatenate((pos_arr,ox_arr,pattern_arr))
        elif N_fing == 38:
            fprint = np.concatenate((pos_arr,ox_arr,dihedral_arr))
        elif N_fing >= 40:
            fprint = np.concatenate((pos_arr,ox_arr,pattern_arr,dihedral_arr))

        geom_mat[i_file] = fprint
        i_file=i_file+1

    Int_lam[N_file,:] = lambda_min
    Int_lam[N_file+1,:] = dlambda
    np.savetxt(file_spec, Int_lam)
    np.savetxt(file_geom, geom_mat)
    return Int_lam, lambda_min, dlambda

