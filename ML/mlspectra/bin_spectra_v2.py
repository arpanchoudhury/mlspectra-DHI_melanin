import pandas as pd
import numpy as np
import mlspectra

def bin_spectra_waveuniform_v2(file_geom,  N_fing, spec_files, dihedral_file):
    #spec_files = mlspectra.read_files(spec_path)
    fo = open(spec_files, 'r')
    filenames = fo.readlines()

    dihedrals = np.loadtxt(dihedral_file)
    dihedrals = np.abs(dihedrals)

    #e_max = 6.1992
    #e_min = 1.5498
    #N_bin=30


    N_file = len(filenames)
    i_file = 0
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
	print(i_file,"done")
    np.savetxt(file_geom, geom_mat)
    return

