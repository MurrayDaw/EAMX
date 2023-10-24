#define parameters for model and call Make_LAMMPS_eamalloy_v2 to write to file
#Daw & Chandross: 8 Dec 2021

from Make_LAMMPS_eamalloy_v2 import Make_LAMMPS_eamalloy_v2
import numpy as np

# define:
#  headers (2)
#  ntypes
#  for each type:
#     beta, r1nn, rho0, phi0, Ec, B, mass, name, rcut
#  for each pair:
#     chi

header0 = '1st line of info'
header1 = '2nd line of info'
headers = [header0,header1]

ntypes = 2
betas = np.ones(ntypes)*1.9
r1nns = np.ones(ntypes)*2.7
rho0s = np.ones(ntypes)*1.0
phi0s = np.ones(ntypes)*0.25
Ecs = np.ones(ntypes)*4.1
Bs = np.ones(ntypes)*1.1
masses = np.ones(ntypes)*30.
names = np.array(['A','B'])
rcuts = np.ones(ntypes)*5.3

chis = np.ones((ntypes,ntypes))

# change lattice constant of second type and scale B, beta accordingly
scale = 1.03
r1nns[1] = r1nns[1]*scale
betas[1] = betas[1]/scale
Bs[1] = Bs[1]/scale**3
rcuts[1] = rcuts[1]*scale

Make_LAMMPS_eamalloy_v2(headers,ntypes,betas,r1nns,rho0s,phi0s,Ecs,Bs,masses,names,rcuts,chis)

print(" perfetto ")
