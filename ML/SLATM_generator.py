import os
import sys
import numpy as np
import qml
from qml.representations  import get_slatm_mbtypes

# Geometry file names and their location
Dir = r"optimized_geometry/"
Ext = r".xyz"

path = os.path.join(os.path.dirname(os.path.realpath(__file__)), Dir)

print(path)
# All molecules
with open("qc_Opt.dat", "r") as fo:
    files = fo.readlines()

mols = []
for xyz_file in files:
    mol = qml.Compound(xyz = path + xyz_file.strip())
    mols.append(mol)

Nmols=len(mols)
print("\n")
print("Loading geometries and properties\n")
print("No. of molecules in the dataset:", Nmols,"\n")

# Many-body types
mbtypes = get_slatm_mbtypes(np.array([mol.nuclear_charges for mol in mols]))
print("Many-body types in the SLATM descriptor:\n",mbtypes,"\n")

# Descriptors
print("Generating slatm descriptor vectors\n")
desc=[]
i=0
for mol in mols:
    mol.generate_slatm(mbtypes, local=False, sigmas=[0.05, 0.05], dgrids=[0.03, 0.03], rcut=5.0, pbc='000', alchemy=False, rpower=6)
    X=np.array(mol.representation)
    desc.append(X)
    i=i+1
    #print(" ...",i," molecules done out of ", Nmols)
desc=np.array(desc)

np.savetxt("SLATM.dat", desc)
