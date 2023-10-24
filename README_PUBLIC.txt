***24 Oct 2023***

This has the EAM-X model in python and also a script to write out a
LAMMPS potential file (eam/alloy format) as well as a simple test.
The python script model_v2.py also does a check on the parameter values
to make sure that the embedding function has positive curvature
everywhere (for positive argument). 

SetParamsWriteFile.py defines parameters for a binary (of types A and
B), calls Make_LAMMPS_eamalloy_v2.py to write the LAMMPS potential
file, based on model_v2.py.

In the example, the two types differ only by lattice constant.

These scripts are all in python3.

Use SetParamsWriteFile.py as a model for how to call Make_LAMMPS_eam_alloya_v2.py.
(That in turn will call model_v2.py.)

example:

python SetParamsWriteFile.py

Then run LAMMPS to evaluate the energy/atom and lattice constant for each type.

lmps < in.example

Compare your results to out.example

...MSDaw (24 Oct 2023)


P.S. See and reference these two seminal papers:
 (1) M. S. Daw & M. E. Chandross, Acta Materialia, v248 a118771 (2023)
   "Simple parameterization of embedded atom method potentials for FCC metals"
    https://doi.org/10.1016/j.actamat.2023.118771
 (2) M. S. Daw & M. E. Chandross, Acta Materialia, v248 a118771 (2023)
   "Simple parameterization of embedded atom method potentials for FCC alloys"
    https://doi.org/10.1016/j.actamat.2023.118772

P.P.S. BUGFIX: rcut problem was fixed (18 Apr 2023) so now Ece, r1nne, Be,
and U3 will all be correct for rcut out to Sqrt[5]r1nne; beyond that,
paramOK will flag.







