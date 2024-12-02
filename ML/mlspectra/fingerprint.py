import numpy as np


"""oh = {'DHI':np.array([0,0]), 
      'MKI':np.array([1,0]),
      'DKI':np.array([1,1])
     }
"""
oh = {'DHI':np.array([1,0,0]), 
      'MKI':np.array([0,1,0]),
      'DKI':np.array([0,0,1])
     }
pos = {'1':np.array([1,0,0,0,0]), 
       '2':np.array([0,1,0,0,0]), 
       '3':np.array([0,0,1,0,0]), 
       '4':np.array([0,0,0,1,0]), 
       '7':np.array([0,0,0,0,1])
      }

def calfingerprint(filename, dihedral_arr_row, discrete):
    mol_name = filename[:15]
    mol_connect = filename[16:27]
    bran = filename[-19:-13]
    ox_list = []
    pos_list = []
    dihedral_list = []
    if bran == 'branch':
        branch = True
    else:
        branch = False

    if branch is True:
        nametoks = mol_name.split('-')
        for i in range(4):      # for tetramer
            ox_list.append(oh[nametoks[i]])


        connect_A = pos[mol_connect[0]]
        pos_list.append(connect_A)

        connect_B = pos[mol_connect[2]] + pos[mol_connect[4]] + pos[mol_connect[8]]
        pos_list.append(connect_B)

        connect_C = pos[mol_connect[6]]
        pos_list.append(connect_C)

        connect_D = pos[mol_connect[-1]]
        pos_list.append(connect_D)

        pattern_arr = np.array([0,1])


    else:
        nametoks = mol_name.split('-')
        for i in range(4):	# for tetramer
            ox_list.append(oh[nametoks[i]])

        connect_A = pos[mol_connect[0]]
        pos_list.append(connect_A)

        connect_B = pos[mol_connect[2]] + pos[mol_connect[4]]
        pos_list.append(connect_B)

        connect_C = pos[mol_connect[6]] + pos[mol_connect[8]]
        pos_list.append(connect_C)

        connect_D = pos[mol_connect[-1]]
        pos_list.append(connect_D)

        pattern_arr = np.array([1,0])
	

    if discrete == 2:
        for j in range(dihedral_arr_row.shape[0]):
            if dihedral_arr_row[j] <= 90.0:
                dihedral_list.append(np.array([1,0]))
            else:
                dihedral_list.append(np.array([0,1]))
    elif discrete == 4:
        for j in range(dihedral_arr_row.shape[0]):
            if dihedral_arr_row[j] <= 45.0:
                dihedral_list.append(np.array([1,0,0,0]))
            elif (45.0 < dihedral_arr_row[j] <= 90.0):
                dihedral_list.append(np.array([0,1,0,0]))
            elif (90.0 < dihedral_arr_row[j] <= 135.0):
                dihedral_list.append(np.array([0,0,1,0]))
            elif (135.0 < dihedral_arr_row[j] <= 180.0):
                dihedral_list.append(np.array([0,0,0,1]))
    elif discrete == 10:
        for j in range(dihedral_arr_row.shape[0]):
            if dihedral_arr_row[j] <= 18.0:
                dihedral_list.append(np.array([1,0,0,0,0,0,0,0,0,0]))
            elif (18.0 < dihedral_arr_row[j] <= 36.0):
                dihedral_list.append(np.array([0,1,0,0,0,0,0,0,0,0]))
            elif (36.0 < dihedral_arr_row[j] <= 54.0):
                dihedral_list.append(np.array([0,0,1,0,0,0,0,0,0,0]))
            elif (54.0 < dihedral_arr_row[j] <= 72.0):
                dihedral_list.append(np.array([0,0,0,1,0,0,0,0,0,0]))
            elif (72.0 < dihedral_arr_row[j] <= 90.0):
                dihedral_list.append(np.array([0,0,0,0,1,0,0,0,0,0]))
            elif (90.0 < dihedral_arr_row[j] <= 108.0):
                dihedral_list.append(np.array([0,0,0,0,0,1,0,0,0,0]))
            elif (108.0 < dihedral_arr_row[j] <= 126.0):
                dihedral_list.append(np.array([0,0,0,0,0,0,1,0,0,0]))
            elif (126.0 < dihedral_arr_row[j] <= 144.0):
                dihedral_list.append(np.array([0,0,0,0,0,0,0,1,0,0]))
            elif (144.0 < dihedral_arr_row[j] <= 162.0):
                dihedral_list.append(np.array([0,0,0,0,0,0,0,0,1,0]))
            elif (162.0 < dihedral_arr_row[j] <= 180.0):
                dihedral_list.append(np.array([0,0,0,0,0,0,0,0,0,1]))
    elif discrete == 0:
        dihedral_list.append(np.array([0]))

    pos_arr = np.array(pos_list).flatten()
    ox_arr = np.array(ox_list).flatten()
    dihedral_arr = np.array(dihedral_list).flatten()

    print(pos_arr)
    print(ox_arr) 
    print(pattern_arr) 
    return pos_arr, ox_arr, pattern_arr, dihedral_arr

