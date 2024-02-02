import numpy as np
import pandas as pd



df1 = pd.read_csv('DHI2.xyz', sep='\s+', header=None, skiprows=2)
atsym_DHI = (df1.iloc[:,0]).to_numpy()
xyz_DHI = (df1.iloc[:,[1,2,3]]).to_numpy()

"""df1 = pd.read_csv('substituents/mki.xyz', sep='\s+', header=None, skiprows=2)
atsym_MKI = (df1.iloc[:,0]).to_numpy()
xyz_MKI = (df1.iloc[:,[1,2,3]]).to_numpy()

df1 = pd.read_csv('substituents/dki.xyz', sep='\s+', header=None, skiprows=2)
atsym_DKI = (df1.iloc[:,0]).to_numpy()
xyz_DKI = (df1.iloc[:,[1,2,3]]).to_numpy()"""


df2 = pd.read_csv('main.xyz', sep='\s+', header=None, skiprows=2)
atsym_main = (df2.iloc[:,0]).to_numpy()
xyz_main = (df2.iloc[:,[1,2,3]]).to_numpy()


def write_mat(total_atoms, atsym, xyz):
    for i in range(total_atoms):
        line = str(atsym[i][0])+' '
        for j in range(3):
            line = line + '%11.7f'%xyz[i,j]+ '  '
        result.write(str(line) + '\n')



def align_to_Zaxis(xyz, inner, outer):
    xyz = xyz - xyz[inner]

    U1 = xyz[outer]
    Vz = np.array([0, 0, 1.0])

    print("U1", U1)
    U1_proj_YZ = np.array([0,U1[1],U1[2]])
    print("U1_proj_YZ", U1_proj_YZ)

    A = np.degrees(np.arccos((np.dot(U1_proj_YZ,Vz))/(np.linalg.norm(U1_proj_YZ) * np.linalg.norm(Vz))))
    if U1[1] > 0:
        Rx = np.array([[1,0,0],[0,np.cos(np.radians(A)),-np.sin(np.radians(A))],[0,np.sin(np.radians(A)),np.cos(np.radians(A))]])
        xyz = np.transpose(np.dot(Rx,np.transpose(xyz)))
    elif U1[1] < 0:
        Rx = np.array([[1,0,0],[0,np.cos(np.radians(A)),np.sin(np.radians(A))],[0,-np.sin(np.radians(A)),np.cos(np.radians(A))]])
        xyz = np.transpose(np.dot(Rx,np.transpose(xyz)))

    U2 = xyz[outer]
    B = np.degrees(np.arccos((np.dot(U2,Vz))/(np.linalg.norm(U2) * np.linalg.norm(Vz))))
    if U2[0] > 0:
        Ry = np.array([[np.cos(np.radians(B)),0,-np.sin(np.radians(B))],[0,1,0],[np.sin(np.radians(B)),0,np.cos(np.radians(B))]]) 
        xyz = np.transpose(np.dot(Ry,np.transpose(xyz)))
    if U2[0] < 0:
        Ry = np.array([[np.cos(np.radians(B)),0,np.sin(np.radians(B))],[0,1,0],[-np.sin(np.radians(B)),0,np.cos(np.radians(B))]]) 
        xyz = np.transpose(np.dot(Ry,np.transpose(xyz)))


    # align to XZ plane
    temp_vec = xyz[3] - xyz[2]
    Vx = np.array([1.0, 0, 0])

    C = np.degrees(np.arccos((np.dot(temp_vec,Vx))/(np.linalg.norm(temp_vec) * np.linalg.norm(Vx))))
    if U1[1] > 0:
        Rz = np.array([[np.cos(np.radians(C)),-np.sin(np.radians(C)),0],[np.sin(np.radians(C)),np.cos(np.radians(C)),0],[0,0,1]])
        xyz = np.transpose(np.dot(Rz,np.transpose(xyz)))
    elif U1[1] < 0:
        Rz = np.array([[np.cos(np.radians(C)),np.sin(np.radians(C)),0],[-np.sin(np.radians(C)),np.cos(np.radians(C)),0],[0,0,1]])
        xyz = np.transpose(np.dot(Rz,np.transpose(xyz)))


    return xyz



# mention the atoms through which new bond is to be formed
inner_atom = 0
outer_atom = 17
xyz_main_new = align_to_Zaxis(xyz_main, inner_atom, outer_atom)

if xyz_main_new[outer_atom,2] < 0.0:
    Rx = np.array([[1,0,0],[0,np.cos(np.radians(180)),-np.sin(np.radians(180))],[0,np.sin(np.radians(180)),np.cos(np.radians(180))]])
    xyz_main_new = np.transpose(np.dot(Rx,np.transpose(xyz_main_new)))


#----------- can be used loop
xyz_main_new = np.delete(xyz_main_new, outer_atom, 0)
atsym_main_new = np.delete(atsym_main, outer_atom, 0)

c = np.concatenate((xyz_main_new,xyz_DHI), axis=0)
d = np.concatenate((atsym_main_new,atsym_DHI), axis=0)

result = open('dimer_'+str(inner_atom)+str(outer_atom)+'.xyz', 'w+')
result.write(str(c.shape[0])+'\n')
result.write('\n')
write_mat(c.shape[0], d, c)
result.close()
