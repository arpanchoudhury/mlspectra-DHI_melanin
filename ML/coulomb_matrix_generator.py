import numpy as np
import math as ma
import pandas as pd


sorted_cm = open('coulombMatrix.dat', 'w+')


Z = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F':9, 'S':16}


def coulomb_matrix(atomic_sym, x, y, z):
    """
    This function takes in xyz coordinates of a molecule to compute the corresponding coulomb Matrix 
    """
    total_no_of_atoms = len(atomic_sym)
    cm = np.zeros((66,66))
    atomic_sym = np.array(atomic_sym)    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    x = x/0.52917721
    y = y/0.52917721
    z = z/0.52917721
    for i in range(total_no_of_atoms):
        for j in range(total_no_of_atoms):
            if i==j:
                cm[i,j] = 0.5*(Z[atomic_sym[i]])**2.4   # Diagonal term described by Potential energy of isolated atom
            else:
                cm[i,j] = (Z[atomic_sym[i]]*Z[atomic_sym[j]])/ma.sqrt(
                          (x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2)   # Pair-wise repulsion 

    #unsorted_cm.write(str(cm) + 2*'\n')
    eigen_vals,v = np.linalg.eig(cm)
    idx = np.argsort(eigen_vals)[::-1]
    sorted_eigen_vals = eigen_vals[idx]
    # writes the unsorted coulomb matrix 
    row_norms = np.linalg.norm(cm, axis=1)
    # calculates the norm of each row of the coulomb matrix
    descending_row_idx = np.argsort(row_norms)[::-1]
    # evaluates the row indices in descending order
    #print(descending_row_idx)
    shuffled_cm = cm[descending_row_idx,:]
    shuffled_cm = shuffled_cm[:,descending_row_idx]
    # sorts the coulomb matrix by shuffling the rows & cols to arrange it in 
    # descending row norms
    #print(np.linalg.norm(shuffled_cm, axis=1))
    lower_trin = shuffled_cm[np.tril_indices(shuffled_cm.shape[0], k = -1)]
    # extracts the lower triangular part along with diag elements of the sorted 
    # coulomb matrix and put them in a 1D array
    for i in range(lower_trin.shape[0]):
        if i < (lower_trin.shape[0]-1):
            sorted_cm.write(str(lower_trin[i])+' ')
        else:
            sorted_cm.write(str(lower_trin[i]))
    sorted_cm.write("\n")
    print(lower_trin.shape[0])
    #lower_trin = lower_trin[:, np.newaxis]
    # converts the 1D array to a column vector
    return lower_trin


# Reading geometry from XYZ file
with open("qc_Opt.dat", "r") as fo:
    files = fo.readlines()

for geom in (files):
    df = pd.read_csv("optimized_geometry/"+geom.strip(), skiprows=2, header=None, sep="\s+")
    x = df.iloc[:,1].to_numpy()
    y = df.iloc[:,2].to_numpy()
    z = df.iloc[:,3].to_numpy()
    atomic_sym = df.iloc[:,0].to_numpy()

    coulomb_matrix(atomic_sym, x, y, z)

sorted_cm.close()

