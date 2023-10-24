#input EAM-X model and parameters, write an LAMMPS eam/alloy file
#Daw & Chandross: 8 Dec 2021
#checks parameters of each element for negative curvature

def Make_LAMMPS_eamalloy_v2(headers,ntypes,betas,r1nns,rho0s,phi0s,Ecs,Bs,masses,names,rcuts,chis):

#    print(" inside Make_LAMMPS_eamalloy_v2 ")
#    print(" headers = ")
#    print(headers)
#    print(" ntypes = ",ntypes)
#    print(" betas = ",betas)
#    print(" r1nns = ",r1nns)
#    print(" rho0s = ",rho0s)
#    print(" phi0s = ",phi0s)
#    print(" Ecs = ",Ecs)
#    print(" Bs = ",Bs)
#    print(" masses = ",masses)
#    print(" names = ",names)
#    print(" rcuts = ",rcuts)
#    print(" chis = ",chis)
    
    from model_v2 import rho, phi, F, paramOK
    import numpy as np
    import math
    import os

    header1 = headers[0]
    header2 = headers[1]

# define r grid
    nr = 10000
    rcutout = ((nr+10.)/nr)*np.max(rcuts)
#    print(" rcuts = ",rcuts)
#    print(" rcut for file = ",rcutout)
    drout = rcutout/nr
    rs = drout*np.array([ir for ir in range(nr)]) 
    rs[0] = rs[1]

# define rho grid
    r1closest = 0.6*np.min(r1nns)
    r2closest = math.sqrt(2.)*r1closest
    r3closest = math.sqrt(3.)*r1closest
    r4closest = math.sqrt(4.)*r1closest
    rhobarmaxes = np.zeros(ntypes)
    for itype in range(ntypes):
        rhobarmaxes[itype] = 12*rho(r1closest,rho0s[itype],betas[itype],r1nns[itype],rcuts[itype]) + \
          6*rho(r2closest,rho0s[itype],betas[itype],r1nns[itype],rcuts[itype]) + \
          24*rho(r3closest,rho0s[itype],betas[itype],r1nns[itype],rcuts[itype]) + \
          12*rho(r4closest,rho0s[itype],betas[itype],r1nns[itype],rcuts[itype])

    rhobarmax = np.max(rhobarmaxes)                        
    nrho = 10000
    drho = rhobarmax/nrho
    rhos = drho*np.array([irho for irho in range(nrho)])

    Z = 1

# define name of file; delete if it already exists
    fileout = "EAMX.eamalloy"
    file_exists = os.path.exists(fileout)
    if file_exists:
        print(fileout," exists already; will delete now ")
        os.remove(fileout)
    
# check that the params for each type are OK
# if not, then quit now (so that no file will exist)
    allOK = 1
    for itype in range(ntypes):
        test = paramOK(phi0s[itype],rho0s[itype],betas[itype],r1nns[itype],Bs[itype],Ecs[itype],rcuts[itype])
        if test == 0:
            print(" paramOK failed for type ",itype)
#        else:
#            print(" paramOK OK for type ",itype)
        allOK = allOK * test
    if allOK == 0:
        print(" somebody failed: exit here w/o starting new file ")
        exit()
#    else:
#        print(" everybody OK ")
        
# OK, now write out the file in LAMMPS eam/alloy format
    print(" starting to write ",fileout)
    stream = open(fileout,"w")    

    header1out = "writing EAM-X functions from model and parameters to LAMMPS eam/alloy \n"
    header2out = header1+"\n"
    header3out = header2+"\n"
    allnames = ""
    for itype in range(ntypes):
        allnames = allnames+" "+names[itype]
    header4out = str(ntypes)+" "+allnames+" \n" # ntypes elements in this file

    stream.write(header1out)
    stream.write(header2out)
    stream.write(header3out)
    stream.write(header4out)

    line5out = " "+str(nrho)+" "+str(drho)+" "+str(nr)+" "+str(drout)+" "+str(rcutout)+" \n"
    stream.write(line5out)

# write rho(r), F(rho) for each element 
    for itype in range(ntypes):
        alat = r1nns[itype]*math.sqrt(2.)
        nextline = str(Z)+" "+str(masses[itype])+" "+str(int(alat*100000)/100000.)+"  FCC \n"
        stream.write(nextline)

# evaluate functions on grids
        Fofrho = F(rhos,phi0s[itype],rho0s[itype],betas[itype],r1nns[itype],Bs[itype],Ecs[itype],rcuts[itype])
        np.savetxt(stream,Fofrho)
        rhoofr = rho(rs,rho0s[itype],betas[itype],r1nns[itype],rcuts[itype])
        np.savetxt(stream,rhoofr)

    for itype in range(ntypes):
        phiiofr = phi(rs,phi0s[itype],betas[itype],r1nns[itype],rcuts[itype])
        for jtype in range(itype+1):
            phijofr = phi(rs,phi0s[jtype],betas[jtype],r1nns[jtype],rcuts[jtype])
            phiijofr = chis[itype,jtype]*(phiiofr + phijofr)/2.
            rphiofr = rs*phiijofr
            np.savetxt(stream,rphiofr)
#            print(" itype, jtype = ",itype,jtype)

    stream.close()
#    print(" perfetto ")        


   



    

