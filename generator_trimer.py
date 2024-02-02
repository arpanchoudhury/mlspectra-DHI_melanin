import os
import numpy as np
import pandas as pd


# atom position dict of substituents
pos_dict_DHI = {0:"2A",1:"3A",4:"1A",5:"7B",10:"4B"}
pos_dict_MKI = {0:"2A", 1:"3A", 10:"4B", 5:"7B"}
pos_dict_DKI = {4:"1A", 0:"2A", 1:"3A", 10:"4B", 5:"7B"}

# main molecules
starting_geoms = "dimer"
resulting_geoms = "trimer"
namestringlist = []

data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), starting_geoms)
all_files = os.listdir(data_dir)




def align_to_Zaxis(xyz, inner, outer, substituent):
    xyz = xyz - xyz[inner]

    U1 = xyz[outer]
    Vz = np.array([0, 0, 1.0])

    U1_proj_YZ = np.array([0,U1[1],U1[2]])

    A = np.degrees(np.arccos((np.dot(U1_proj_YZ,Vz))/(np.linalg.norm(U1_proj_YZ) * np.linalg.norm(Vz))))
    if U1[1] > 0:
        Rx = np.array([[1,0,0],
                       [0,np.cos(np.radians(A)),-np.sin(np.radians(A))],
                       [0,np.sin(np.radians(A)),np.cos(np.radians(A))]])
        xyz = np.transpose(np.dot(Rx,np.transpose(xyz)))
    elif U1[1] < 0:
        Rx = np.array([[1,0,0],
                       [0,np.cos(np.radians(A)),np.sin(np.radians(A))],
                       [0,-np.sin(np.radians(A)),np.cos(np.radians(A))]])
        xyz = np.transpose(np.dot(Rx,np.transpose(xyz)))

    U2 = xyz[outer]
    B = np.degrees(np.arccos((np.dot(U2,Vz))/(np.linalg.norm(U2) * np.linalg.norm(Vz))))
    if U2[0] > 0:
        Ry = np.array([[np.cos(np.radians(B)),0,-np.sin(np.radians(B))],
                       [0,1,0],
                       [np.sin(np.radians(B)),0,np.cos(np.radians(B))]]) 
        xyz = np.transpose(np.dot(Ry,np.transpose(xyz)))
    if U2[0] < 0:
        Ry = np.array([[np.cos(np.radians(B)),0,np.sin(np.radians(B))],
                       [0,1,0],
                       [-np.sin(np.radians(B)),0,np.cos(np.radians(B))]]) 
        xyz = np.transpose(np.dot(Ry,np.transpose(xyz)))


    # align to XZ plane
    dist_list_idx = []
    for i in range(xyz.shape[0]):
        dist = np.linalg.norm(xyz[i] - xyz[inner]) 
        if (dist <= 1.5) and (i not in [inner,outer]):
            dist_list_idx.append(i)
    
    print("dist_list_idx[0]",dist_list_idx[0])
    print("dist_list_idx[1]",dist_list_idx[1])
    temp_vec1 = xyz[dist_list_idx[0]] - xyz[inner]
    temp_vec2 = xyz[dist_list_idx[1]] - xyz[inner]
    normal = np.cross(temp_vec1, temp_vec2) 
    Vx = np.array([1.0, 0, 0])

    #d = normal[0]*xyz[2,0] + normal[1]*xyz[2,1] + normal[2]*xyz[2,2]
    C = np.degrees(np.arccos((np.dot(normal,Vx))/(np.linalg.norm(normal) * np.linalg.norm(Vx))))
    C = 90.0 - C


    """for i in range(xyz.shape[0]):
        if xyz[i,0] > 0.0 and i not in [inner,outer]:
            tmp_idx = i
            break"""

    

    for i in dist_list_idx:
        if xyz[i,0] > 0.0:
            tmp_idx = i
            break  
        elif xyz[i,0] == 0.0:
            tmp_idx = None
            break

    """if substituent:      
        print(tmp_idx)
        print("angle =", C)"""
    
    """if (tmp_idx != None) and (xyz[tmp_idx,1] < 0):
        Rz = np.array([[np.cos(np.radians(C)),-np.sin(np.radians(C)),0],
                       [np.sin(np.radians(C)),np.cos(np.radians(C)),0],
                       [0,0,1]])
        xyz = np.transpose(np.dot(Rz,np.transpose(xyz)))
    elif (tmp_idx != None) and (xyz[tmp_idx,1] > 0):
        Rz = np.array([[np.cos(np.radians(C)),np.sin(np.radians(C)),0],
                       [-np.sin(np.radians(C)),np.cos(np.radians(C)),0],
                       [0,0,1]])
        xyz = np.transpose(np.dot(Rz,np.transpose(xyz)))
    elif tmp_idx == None:
        Rz = np.array([[np.cos(np.radians(90)),np.sin(np.radians(90)),0],
                       [-np.sin(np.radians(90)),np.cos(np.radians(90)),0],
                       [0,0,1]])
        xyz = np.transpose(np.dot(Rz,np.transpose(xyz)))"""


    Rz = np.array([[np.cos(np.radians(C)),-np.sin(np.radians(C)),0],
                       [np.sin(np.radians(C)),np.cos(np.radians(C)),0],
                       [0,0,1]])
    xyz = np.transpose(np.dot(Rz,np.transpose(xyz)))    

    if substituent:
        Rx = np.array([[1,0,0],
                       [0,np.cos(np.radians(180)),-np.sin(np.radians(180))],
                       [0,np.sin(np.radians(180)),np.cos(np.radians(180))]])
        xyz = np.transpose(np.dot(Rx,np.transpose(xyz)))
        xyz[:,2] = xyz[:,2] + 1.47

    return xyz


# mention the atoms through which new bond is to be formed
# ----- main molecule -----
for file in all_files:
#print(all_files[2])
#for file in [all_files[1]]:
    print(file)
    data_file = os.path.join(data_dir, file)
    df = pd.read_csv(data_file, sep='\s+', header=None, skiprows=1)
    avail_pos = (df.iloc[0]).to_list()

    tmp_len = int((len(avail_pos)-5)/3) #
    inner_atom_main_list = [i[:-1] for i in avail_pos[:tmp_len]]
    inner_ring_main_list = [i[-1] for i in avail_pos[:tmp_len]]
    outer_atom_main_list = [i for i in avail_pos[tmp_len:int(2*tmp_len)]]
    inner_pos_main_list = [i for i in avail_pos[int(2*tmp_len):int(3*tmp_len)]]
    n_subs_per_ring = [i for i in avail_pos[int(3*tmp_len):-1]]
    symmetry = avail_pos[-1]

                

    #elif symmetry == "nosym":
    #if file == "DHI-DHI_2-3.xyz":
        #print(inner_atom_main_list)
    inner_atom_main_list = list(map(int, inner_atom_main_list))
    outer_atom_main_list = list(map(int, outer_atom_main_list))
    inner_pos_main_list = list(map(int, inner_pos_main_list))
    n_subs_per_ring = list(map(int, n_subs_per_ring))

    frag_atom = []
    frag_pos = []
    m = 0
    n = n_subs_per_ring[0]
    for i in range(len(n_subs_per_ring)):
        if n_subs_per_ring[i] != 0:
            tmp1 = [j for j in inner_atom_main_list[m:n]]
            tmp2 = [j for j in inner_pos_main_list[m:n]]
            frag_atom.append(tmp1)
            frag_pos.append(tmp2)
        elif n_subs_per_ring[i] == 0:
            frag_atom.append([])
            frag_pos.append([])


        if i == len(n_subs_per_ring)-1:
            break

        m = m + n_subs_per_ring[i]
        n = n + n_subs_per_ring[i+1]


    #if file == "DHI-DHI_2-3.xyz":
        #print(frag_atom)
    atsym_main = (df.iloc[1:,0]).to_numpy()
    xyz_main = (df.iloc[1:,[1,2,3]]).to_numpy(dtype=np.float32)
    ## -----

    for idx_main in range(len(inner_atom_main_list)):

        inner_atom_main = inner_atom_main_list[idx_main]
        outer_atom_main = outer_atom_main_list[idx_main]
        print("inner_atom_main",inner_atom_main)
        ring_type = inner_ring_main_list[idx_main]


        substituent = False
        xyz_main_copy = xyz_main.copy()
        atsym_main_copy = atsym_main.copy()

        xyz_main_new = align_to_Zaxis(xyz_main_copy, inner_atom_main, outer_atom_main, substituent)
        if xyz_main_new[outer_atom_main,2] < 0.0:
            Rx = np.array([[1,0,0],
                           [0,np.cos(np.radians(180)),-np.sin(np.radians(180))],
                           [0,np.sin(np.radians(180)),np.cos(np.radians(180))]])
            xyz_main_new = np.transpose(np.dot(Rx,np.transpose(xyz_main_new)))

        

        xyz_main_new = np.delete(xyz_main_new, outer_atom_main, 0)
        atsym_main_new = np.delete(atsym_main_copy, outer_atom_main, 0)
        # ----- X -----


        # loop over differnt type of substituents
        for subs in ["DHI", "DKI", "MKI"]:
            
            substituent = True

            df1 = pd.read_csv("substituents/"+subs+".xyz", sep='\s+', header=None, skiprows=2)
            atsym_subs = (df1.iloc[:,0]).to_numpy()
            xyz_subs = (df1.iloc[:,[1,2,3]]).to_numpy()

            # load index file to store the available indices for oligomerization
            index1 = np.loadtxt("index1_"+subs+".txt", dtype="int") # for 5-mem ring
            index2 = np.loadtxt("index2_"+subs+".txt", dtype="int") # for 6-mem ring

            # 
            """if ring_type == "B":
                tmp_inner_list = list(index1[0])
                tmp_outer_list = list(index1[1])
            else:
                tmp_inner_list = list(index1[0]) + list(index2[0])
                tmp_outer_list = list(index1[1]) + list(index2[1])"""
            tmp_inner_list = list(index1[0]) + list(index2[0])
            tmp_outer_list = list(index1[1]) + list(index2[1])


            #print("inner_atom_main_list[idx_main]", inner_atom_main_list[idx_main])
            for i in range(len(frag_atom)):
                if inner_atom_main_list[idx_main] in frag_atom[i]:
                    tmp_list_atom = frag_atom[i].copy()
                    tmp_list_pos = frag_pos[i].copy()
                    substituted_ring = i
            #tmp_list.remove(inner_pos_main_list[idx_main])
            #print("tmp_list_atom",tmp_list_atom)
            #print("tmp_list_pos",tmp_list_pos)
            """if inner_pos_main_list[idx_main] == 1:
                tmp_pos_main = [3,4]
            elif inner_pos_main_list[idx_main] == 2:
                tmp_pos_main = [4,7]   
            elif inner_pos_main_list[idx_main] == 3 and subs != "MKI":
                tmp_pos_main = [1,7]
            elif inner_pos_main_list[idx_main] == 3 and subs == "MKI":
                tmp_pos_main = [7]
            elif inner_pos_main_list[idx_main] == 4 and subs != "MKI":
                tmp_pos_main = [2,1,7]  
            elif inner_pos_main_list[idx_main] == 4 and subs == "MKI":
                tmp_pos_main = [2,7]
            elif inner_pos_main_list[idx_main] == 7:
                tmp_pos_main = [2,3,4]"""       

            if inner_pos_main_list[idx_main] == 1:
                tmp_pos_main = [3,4]
            elif inner_pos_main_list[idx_main] == 2:
                tmp_pos_main = [4,7]   
            elif inner_pos_main_list[idx_main] == 3:
                tmp_pos_main = [1,7]
            elif inner_pos_main_list[idx_main] == 4:
                tmp_pos_main = [2,1,7]  
            elif inner_pos_main_list[idx_main] == 7:
                tmp_pos_main = [2,3,4]
            
            #print("inner_pos_main_list",inner_pos_main_list)
            #print("tmp_pos_main",tmp_pos_main)
            tmp_i = []
            for i in range(len(tmp_list_pos)):
                if tmp_list_pos[i] != inner_pos_main_list[idx_main] and tmp_list_pos[i] in tmp_pos_main:
                    #tmp_i.append(inner_pos_main_list.index(tmp_list_pos[i]))
                    tmp_i.append(inner_atom_main_list.index(tmp_list_atom[i]))



            remaining_idx_list = []
            for i in range(len(frag_atom)):
                if tmp_list_atom != frag_atom[i]:
                    for j in range(len(frag_atom[i])):
                        remaining_idx_list.append(inner_atom_main_list.index(frag_atom[i][j]))

            #print("remaining_idx_list",remaining_idx_list)
            remaining_idx_list = remaining_idx_list + tmp_i
            #remaining_idx_list.sort()
            #print("inner_atom_main_list",inner_atom_main_list)
            #print("frag_atom",frag_atom)
            #print("frag_pos",frag_pos)
            #print("remaining_idx_list",remaining_idx_list)
            remaining_atom_main_inner = np.array(inner_atom_main_list)[remaining_idx_list]
            #print("remaining_atom_main_inner",remaining_atom_main_inner)
            remaining_pos_main_inner = np.array(inner_pos_main_list)[remaining_idx_list]
            remaining_ring_main_inner = np.array(inner_ring_main_list)[remaining_idx_list]
            remaining_atom_main_outer = np.array(outer_atom_main_list)[remaining_idx_list]
            #if file == "DHI-DHI_2-2.xyz":
                #print(remaining_atom_main_inner)
            for k in range(len(remaining_atom_main_outer)):
                    if remaining_atom_main_outer[k] > outer_atom_main:
                        remaining_atom_main_outer[k] = remaining_atom_main_outer[k] - 1
                    if remaining_atom_main_inner[k] > outer_atom_main:
                        remaining_atom_main_inner[k] = remaining_atom_main_inner[k] - 1

            #if file == "DHI-DHI_2-2.xyz":
                #print(remaining_atom_main_inner)
            # loop over all possible indices available to oligomerize, of the given substituent
            for idx in range(len(tmp_inner_list)):
                xyz_subs_new = align_to_Zaxis(xyz_subs.copy(), tmp_inner_list[idx], tmp_outer_list[idx], substituent)

                xyz_subs_new = np.delete(xyz_subs_new, tmp_outer_list[idx], 0)
                atsym_subs_new = np.delete(atsym_subs.copy(), tmp_outer_list[idx], 0)


                # writing the resulting geometries in XYZ files
                exec("pos_dict = pos_dict_" + subs)
                if resulting_geoms == "dimer":
                    filestring1 = file[:3]
                    filestring2 = ''
                elif resulting_geoms == "trimer":
                    filestring1 = file[:7]
                    filestring2 = file[8:-4]
                    
                elif resulting_geoms == "tetramer":
                    filestring1 = file[:11]
                    filestring2 = file[12:-4]

                b = [int(i) for i in filestring2.split("-") if i.isdigit()]
                for i in range(len(b)):
                    if b[i] > outer_atom_main:
                        b[i] = b[i] - 1
                filestring2 = '-'.join(str(i) for i in b)


                # write the remaining positions of subs (also avoiding steric hindrence)
                if tmp_inner_list[idx] == 0:
                    remaining_atom_subs_inner = [5,10]
                elif tmp_inner_list[idx] == 1 and subs != "MKI":
                    remaining_atom_subs_inner = [4,5]
                elif tmp_inner_list[idx] == 1 and subs == "MKI":
                    remaining_atom_subs_inner = [5]
                elif tmp_inner_list[idx] == 4:
                    remaining_atom_subs_inner = [1,10]
                elif tmp_inner_list[idx] == 5:
                    remaining_atom_subs_inner = [0,1,10]
                elif tmp_inner_list[idx] == 10 and subs != "MKI":
                    remaining_atom_subs_inner = [0,4,5]
                elif tmp_inner_list[idx] == 10 and subs == "MKI":
                    remaining_atom_subs_inner = [0,5]


                remaining_ring_subs_inner = [pos_dict[i][1] for i in remaining_atom_subs_inner]
                remaining_pos_subs_inner = [pos_dict[i][0] for i in remaining_atom_subs_inner]
                list0 = list(index1[0]) + list(index2[0])
                list1 = list(index1[1]) + list(index2[1])
                remaining_atom_subs_outer_idx = [list0.index(i) for i in remaining_atom_subs_inner]
                remaining_atom_subs_outer = [list1[i] for i in remaining_atom_subs_outer_idx]       
                remaining_atom_subs_inner = list(s+xyz_main_new.shape[0] for s in remaining_atom_subs_inner)

                for k in range(len(remaining_atom_subs_outer)):
                    if remaining_atom_subs_outer[k] > tmp_outer_list[idx]:
                        remaining_atom_subs_outer[k] = remaining_atom_subs_outer[k] - 1
                remaining_atom_subs_outer = list(s+xyz_main_new.shape[0] for s in remaining_atom_subs_outer)

                remaining_atom_main_inner = np.array(remaining_atom_main_inner)
                remaining_atom_main_outer = np.array(remaining_atom_main_outer)
                remaining_pos_main_inner = np.array(remaining_pos_main_inner)

                #print("remaining_atom_main_inner",remaining_atom_main_inner)
                remaining_atom_main_inner_idx = np.argsort(remaining_atom_main_inner)
                #print("inner_atom_main",inner_atom_main)
                #print("remaining_atom_main_inner",remaining_atom_main_inner)
                remaining_atom_main_inner = remaining_atom_main_inner[remaining_atom_main_inner_idx].tolist()
                remaining_atom_main_outer = remaining_atom_main_outer[remaining_atom_main_inner_idx].tolist()
                remaining_pos_main_inner = remaining_pos_main_inner[remaining_atom_main_inner_idx].tolist()

                comment_line1 = ' '.join(str(remaining_atom_main_inner[i])+str(remaining_ring_main_inner[i]) for i in range(len(remaining_atom_main_inner)))
                comment_line2 = ' '.join(str(i) for i in remaining_atom_main_outer)
                comment_line3 = ' '.join(str(i) for i in remaining_pos_main_inner)

                comment_line4 = ' '.join(str(remaining_atom_subs_inner[i])+str(remaining_ring_subs_inner[i]) for i in range(len(remaining_ring_subs_inner)))
                comment_line5 = ' '.join(str(i) for i in remaining_atom_subs_outer)
                comment_line6 = ' '.join(str(i) for i in remaining_pos_subs_inner)


                #print(comment_line1)
                n_subs_per_ring_copy = n_subs_per_ring.copy()
                """li = []
                for i in range(len(n_subs_per_ring_copy)):
                    if n_subs_per_ring_copy[i] != 0:
                        li.append(n_subs_per_ring_copy[i])
                        
                for i in range(3,-1,-1):
                    if n_subs_per_ring_copy[i] != 0:
                        last_nonzero = i
                        break

                for i in range(len(n_subs_per_ring)):
                    if substituted_ring == i:
                        n_subs_per_ring_copy[i] = len(tmp_i)

                n_subs_per_ring_copy[last_nonzero+1] = len(remaining_atom_subs_inner)

                comment_line7 = ' '.join(str(i) for i in n_subs_per_ring_copy)"""



                if 2 not in [inner_pos_main_list[idx_main], int(pos_dict[tmp_inner_list[idx]][0])]:
                    Rz = np.array([[np.cos(np.radians(40)),-np.sin(np.radians(40)),0],
                                   [np.sin(np.radians(40)),np.cos(np.radians(40)),0],
                                   [0,0,1]])
                    xyz_subs_new = np.transpose(np.dot(Rz,np.transpose(xyz_subs_new)))

                xyz_final = np.concatenate((xyz_main_new,xyz_subs_new), axis=0)
                atsym_final = np.concatenate((atsym_main_new,atsym_subs_new), axis=0)

                #namestring = [filestring1[:3], filestring1[:3], subs, int(pos_dict[tmp_inner_list[idx]][0]), inner_pos_main_list[idx_main]]

                
                comment_line8 = 'nosym'
                    
                """comment_line = (comment_line1 + ' ' + comment_line4 + ' ' + comment_line2 
                                + ' ' + comment_line5 +' '+ comment_line3 +' '+ comment_line6 
                                +' '+ comment_line7+' '+ comment_line8)"""


                #if namestring not in namestringlist:
                """print(filestring1)
                print(filestring1[:3])
                print(filestring1[4:7])
                print(filestring2)"""

                if symmetry == "nosym":
                    if substituted_ring == 0:
                        namestring = [subs, filestring1[:3], filestring1[4:7], int(pos_dict[tmp_inner_list[idx]][0]), inner_pos_main_list[idx_main],
                                      int(filestring2[0]), int(filestring2[-1])]
                        
                        if namestring not in namestringlist:
                            #if namestring == ["DHI","DHI","DHI",3,4,2,4]:
                                #print(namestring)
                            result = open(resulting_geoms + "/" + subs +"-"+ filestring1 + "_" 
                                + str(pos_dict[tmp_inner_list[idx]][0]) + "-" 
                                + str(inner_pos_main_list[idx_main]) +"-"+ filestring2 +'.xyz', 'w+')
                            

                            print(subs +"-"+ filestring1 + "_" 
                                + str(pos_dict[tmp_inner_list[idx]][0]) + "-" 
                                + str(inner_pos_main_list[idx_main]) +"-"+ filestring2 +'.xyz')
                            result.write(str(xyz_final.shape[0])+'\n')
                            """if file == "DHI-DHI_2-2.xyz":
                                print(file) 
                                print(comment_line4)
                                print(comment_line1)
                                print(remaining_atom_main_inner)"""

                            n_subs_per_ring_copy = [len(remaining_atom_subs_inner), len(tmp_i), n_subs_per_ring_copy[1], 0]   
                            comment_line7 = ' '.join(str(i) for i in n_subs_per_ring_copy) 
                            comment_line = (comment_line4 + ' ' + comment_line1 + ' ' + comment_line5 
                                + ' ' + comment_line2 +' '+ comment_line6 +' '+ comment_line3 
                                +' '+ comment_line7+' '+ comment_line8)
                            result.write(comment_line+"\n")
                
                            for row in range(xyz_final.shape[0]):
                                line = str(atsym_final[row][0])+'  '
                                for col in range(3):
                                    line = line + '%11.7f'%xyz_final[row,col]+ '  '
                                result.write(str(line) + '\n')

                            result.close()
                        namestringlist.append([subs, filestring1[:3], filestring1[4:7], int(pos_dict[tmp_inner_list[idx]][0]), inner_pos_main_list[idx_main],
                                                    int(filestring2[0]), int(filestring2[-1])])
                        namestringlist.append([filestring1[4:7], filestring1[:3], subs, int(filestring2[-1]), int(filestring2[0]), inner_pos_main_list[idx_main], 
                                                    int(pos_dict[tmp_inner_list[idx]][0])])
                            
                    if substituted_ring == 1:
                        namestring = [filestring1[:3], filestring1[4:7], subs, int(filestring2[0]), int(filestring2[-1]),
                                      inner_pos_main_list[idx_main], int(pos_dict[tmp_inner_list[idx]][0])]
                        
                        if namestring not in namestringlist:
                            #if namestring == ["DHI","DHI","DHI",3,4,2,4]:
                                #print(namestring)
                            result = open(resulting_geoms + "/" + filestring1 +"-"+ subs + "_" + filestring2 +"-"
                                    + str(inner_pos_main_list[idx_main]) + "-"
                                    + str(pos_dict[tmp_inner_list[idx]][0])  +'.xyz', 'w+')
                            
                            print(filestring1 +"-"+ subs + "_" + filestring2 +"-"
                                    + str(inner_pos_main_list[idx_main]) + "-"
                                    + str(pos_dict[tmp_inner_list[idx]][0])  +'.xyz')
                            result.write(str(xyz_final.shape[0])+'\n')

                            n_subs_per_ring_copy = [n_subs_per_ring_copy[0], len(tmp_i), len(remaining_atom_subs_inner),  0]   
                            comment_line7 = ' '.join(str(i) for i in n_subs_per_ring_copy)
                            comment_line = (comment_line1 + ' ' + comment_line4 + ' ' + comment_line2 
                                + ' ' + comment_line5 +' '+ comment_line3 +' '+ comment_line6 
                                +' '+ comment_line7+' '+ comment_line8)
                            result.write(comment_line+"\n")
                
                            for row in range(xyz_final.shape[0]):
                                line = str(atsym_final[row][0])+'  '
                                for col in range(3):
                                    line = line + '%11.7f'%xyz_final[row,col]+ '  '
                                result.write(str(line) + '\n')

                            result.close()
                        namestringlist.append([filestring1[:3], filestring1[4:7], subs, int(filestring2[0]), int(filestring2[-1]),
                                                inner_pos_main_list[idx_main], int(pos_dict[tmp_inner_list[idx]][0])])
                        namestringlist.append([subs, filestring1[4:7], filestring1[:3], int(pos_dict[tmp_inner_list[idx]][0]), inner_pos_main_list[idx_main],
                                                int(filestring2[-1]), int(filestring2[0])])
                            
                elif symmetry == "sym":
                    """namestring = [filestring1[:3], filestring1[4:7], subs, int(filestring2[0]), int(filestring2[-1]),
                                  inner_pos_main_list[idx_main], int(pos_dict[tmp_inner_list[idx]][0])]
                    

                    
                    if file == "DHI-DHI_2-2.xyz":
                                print(file) 
                                print(comment_line4)
                                print(comment_line1)
                                print(remaining_atom_main_inner)
                    
                    if namestring not in namestringlist:
                        #if namestring == ["DHI","DHI","DHI",3,4,2,4]:
                            #print(namestring)
                        result = open(resulting_geoms + "/" + filestring1 +"-"+ subs + "_" + filestring2 +"-"
                                    + str(inner_pos_main_list[idx_main]) + "-"
                                    + str(pos_dict[tmp_inner_list[idx]][0]) +'.xyz', 'w+')

                        result.write(str(xyz_final.shape[0])+'\n')

                        n_subs_per_ring_copy = [n_subs_per_ring_copy[0], len(tmp_i), len(remaining_atom_subs_inner),  0]   
                        comment_line7 = ' '.join(str(i) for i in n_subs_per_ring_copy)
                        comment_line = (comment_line1 + ' ' + comment_line4 + ' ' + comment_line2 
                                + ' ' + comment_line5 +' '+ comment_line3 +' '+ comment_line6 
                                +' '+ comment_line7+' '+ comment_line8)
                        result.write(comment_line+"\n")
                
                        for row in range(xyz_final.shape[0]):
                            line = str(atsym_final[row][0])+'  '
                            for col in range(3):
                                line = line + '%11.7f'%xyz_final[row,col]+ '  '
                            result.write(str(line) + '\n')

                        result.close()                        
                        namestringlist.append([filestring1[:3], filestring1[4:7], subs, int(filestring2[0]), int(filestring2[-1]),
                                            inner_pos_main_list[idx_main], int(pos_dict[tmp_inner_list[idx]][0])])
                        namestringlist.append([subs, filestring1[4:7], filestring1[:3], int(pos_dict[tmp_inner_list[idx]][0]),
                                              inner_pos_main_list[idx_main], int(filestring2[-1]), int(filestring2[0])])"""
                    
                    if substituted_ring == 0:
                        namestring = [subs, filestring1[:3], filestring1[4:7], int(pos_dict[tmp_inner_list[idx]][0]), inner_pos_main_list[idx_main],
                                      int(filestring2[0]), int(filestring2[-1])]
                        
                        if namestring not in namestringlist:
                            #if namestring == ["DHI","DHI","DHI",3,4,2,4]:
                                #print(namestring)
                            result = open(resulting_geoms + "/" + subs +"-"+ filestring1 + "_" 
                                + str(pos_dict[tmp_inner_list[idx]][0]) + "-" 
                                + str(inner_pos_main_list[idx_main]) +"-"+ filestring2 +'.xyz', 'w+')
                            
                            result.write(str(xyz_final.shape[0])+'\n')
                            #if file == "DHI-DHI_2-2.xyz":
                                #print(file) 
                                #print(comment_line4)
                                #print(comment_line1)
                                #print(remaining_atom_main_inner)

                            n_subs_per_ring_copy = [len(remaining_atom_subs_inner), len(tmp_i), n_subs_per_ring_copy[1], 0]   
                            comment_line7 = ' '.join(str(i) for i in n_subs_per_ring_copy) 
                            comment_line = (comment_line4 + ' ' + comment_line1 + ' ' + comment_line5 
                                + ' ' + comment_line2 +' '+ comment_line6 +' '+ comment_line3 
                                +' '+ comment_line7+' '+ comment_line8)
                            result.write(comment_line+"\n")
                
                            for row in range(xyz_final.shape[0]):
                                line = str(atsym_final[row][0])+'  '
                                for col in range(3):
                                    line = line + '%11.7f'%xyz_final[row,col]+ '  '
                                result.write(str(line) + '\n')

                            result.close()
                        namestringlist.append([subs, filestring1[:3], filestring1[4:7], int(pos_dict[tmp_inner_list[idx]][0]), inner_pos_main_list[idx_main],
                                                    int(filestring2[0]), int(filestring2[-1])])
                        namestringlist.append([filestring1[4:7], filestring1[:3], subs, int(filestring2[-1]), int(filestring2[0]), inner_pos_main_list[idx_main], 
                                                    int(pos_dict[tmp_inner_list[idx]][0])])
                            
                    if substituted_ring == 1:
                        namestring = [filestring1[:3], filestring1[4:7], subs, int(filestring2[0]), int(filestring2[-1]),
                                      inner_pos_main_list[idx_main], int(pos_dict[tmp_inner_list[idx]][0])]
                        
                        if namestring not in namestringlist:
                            #if namestring == ["DHI","DHI","DHI",3,4,2,4]:
                                #print(namestring)
                            result = open(resulting_geoms + "/" + filestring1 +"-"+ subs + "_" + filestring2 +"-"
                                    + str(inner_pos_main_list[idx_main]) + "-"
                                    + str(pos_dict[tmp_inner_list[idx]][0])  +'.xyz', 'w+')
                            
                            result.write(str(xyz_final.shape[0])+'\n')

                            n_subs_per_ring_copy = [n_subs_per_ring_copy[0], len(tmp_i), len(remaining_atom_subs_inner),  0]   
                            comment_line7 = ' '.join(str(i) for i in n_subs_per_ring_copy)
                            comment_line = (comment_line1 + ' ' + comment_line4 + ' ' + comment_line2 
                                + ' ' + comment_line5 +' '+ comment_line3 +' '+ comment_line6 
                                +' '+ comment_line7+' '+ comment_line8)
                            result.write(comment_line+"\n")
                
                            for row in range(xyz_final.shape[0]):
                                line = str(atsym_final[row][0])+'  '
                                for col in range(3):
                                    line = line + '%11.7f'%xyz_final[row,col]+ '  '
                                result.write(str(line) + '\n')

                            result.close()
                        namestringlist.append([filestring1[:3], filestring1[4:7], subs, int(filestring2[0]), int(filestring2[-1]),
                                                inner_pos_main_list[idx_main], int(pos_dict[tmp_inner_list[idx]][0])])
                        namestringlist.append([subs, filestring1[4:7], filestring1[:3], int(pos_dict[tmp_inner_list[idx]][0]), inner_pos_main_list[idx_main],
                                                int(filestring2[-1]), int(filestring2[0])])

                
                
                

