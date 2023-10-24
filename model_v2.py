# EAMX model with param checker (v2)
# written fully in python (not translated from mathematica)
# Daw and Chandross: 21 Apr 2023
# see and reference these two seminal papers:
# (1) M. S. Daw & M. E. Chandross, Acta Materialia, v248 a118771 (2023)
#   "Simple parameterization of embedded atom method potentials for FCC metals"
#    https://doi.org/10.1016/j.actamat.2023.118771
# (2) M. S. Daw & M. E. Chandross, Acta Materialia, v248 a118771 (2023)
#   "Simple parameterization of embedded atom method potentials for FCC alloys"
#    https://doi.org/10.1016/j.actamat.2023.118772
#
# version 2
# full revision to simply python
# and to fix bug in version 1 (rcut was incorrectly implemented)
#


import numpy as np
import math


# FCC
Zs = np.array([12,6,24,12])
zetas = np.sqrt(np.array([1,2,3,4]))
nshellmax = Zs.size


def fz(z):
    return np.exp(-z)-1+z

def fzp(z):
    return -np.exp(-z)+1

def fzpp(z):
    return np.exp(-z)

def fzppp(z):
    return -np.exp(-z)

# rho(r) and derivatives
# ignore derivatives of cutoff

def cutoff(r,rcut):
    return np.heaviside(rcut-r,0.5)

def rho(r,rho0,beta,r1nne,rcut):
    z = beta*(r-rcut)
    z1 = beta*(r1nne-rcut)
    return rho0*(fz(z)/fz(z1))*cutoff(r,rcut)

def rhop(r,rho0,beta,r1nne,rcut):
    z = beta*(r-rcut)
    z1 = beta*(r1nne-rcut)
    return rho0*(beta*fzp(z)/fz(z1))*cutoff(r,rcut)

def rhopp(r,rho0,beta,r1nne,rcut):
    z = beta*(r-rcut)
    z1 = beta*(r1nne-rcut)
    return rho0*(beta**2*fzpp(z)/fz(z1))*cutoff(r,rcut)

def rhoppp(r,rho0,beta,r1nne,rcut):
    z = beta*(r-rcut)
    z1 = beta*(r1nne-rcut)
    return rho0*(beta**3*fzppp(z)/fz(z1))*cutoff(r,rcut)

# define rhobar and derivatives w.r.t. r1nn 
def rhobar(r1nn,rho0,beta,r1nne,rcut):
    rs = zetas*r1nn
#    rhos = np.array([rho(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhos = rho(rs,rho0,beta,r1nne,rcut)
    return np.dot(Zs,rhos)

def rhobarp(r1nn,rho0,beta,r1nne,rcut):
    rs = zetas*r1nn
#    rhops = np.array([rhop(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhops = rhop(rs,rho0,beta,r1nne,rcut)
    return np.dot(Zs,rhops*zetas)

def rhobarpp(r1nn,rho0,beta,r1nne,rcut):
    rs = zetas*r1nn
#    rhopps = np.array([rhopp(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhopps = rhopp(rs,rho0,beta,r1nne,rcut)
    return np.dot(Zs,rhopps*zetas**2)

def rhobarppp(r1nn,rho0,beta,r1nne,rcut):
    rs = zetas*r1nn
#    rhoppps = np.array([rhoppp(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhoppps = rhoppp(rs,rho0,beta,r1nne,rcut)
    return np.dot(Zs,rhoppps*zetas**3)

# define rhobare and derivatives (rhobar and derivatives evaluated at equilibrium) 
def rhobare(rho0,beta,r1nne,rcut):
    rs = zetas*r1nne
#    rhos = np.array([rho(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhos = rho(rs,rho0,beta,r1nne,rcut) 
    return np.dot(Zs,rhos)

def rhobarpe(rho0,beta,r1nne,rcut):
    rs = zetas*r1nne
#    rhops = np.array([rhop(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhops = rhop(rs,rho0,beta,r1nne,rcut) 
    return np.dot(Zs,rhops*zetas)

def rhobarppe(rho0,beta,r1nne,rcut):
    rs = zetas*r1nne
#    rhopps = np.array([rhopp(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhopps = rhopp(rs,rho0,beta,r1nne,rcut)
    return np.dot(Zs,rhopps*zetas**2)

def rhobarpppe(rho0,beta,r1nne,rcut):
    rs = zetas*r1nne
#    rhoppps = np.array([rhoppp(rs[irs],rho0,beta,r1nne,rcut) for irs in range(nshellmax)])
    rhoppps = rhoppp(rs,rho0,beta,r1nne,rcut) 
    return np.dot(Zs,rhoppps*zetas**3)

# phi(r) and derivatives
# ignore derivatives of cutoff

def phi(r,phi0,beta,r1nne,rcut):
    gamma = 2*beta
    z = gamma*(r-rcut)
    z1 = gamma*(r1nne-rcut)
    calc = phi0*(fz(z)/fz(z1))*cutoff(r,rcut)
    return calc

def phip(r,phi0,beta,r1nne,rcut):
    gamma = 2*beta
    z = gamma*(r-rcut)
    z1 = gamma*(r1nne-rcut)
    return phi0*(gamma*fzp(z)/fz(z1))*cutoff(r,rcut)

def phipp(r,phi0,beta,r1nne,rcut):
    gamma = 2*beta
    z = gamma*(r-rcut)
    z1 = gamma*(r1nne-rcut)
    return phi0*(gamma**2*fzpp(z)/fz(z1))*cutoff(r,rcut)

def phippp(r,phi0,beta,r1nne,rcut):
    gamma = 2*beta
    z = gamma*(r-rcut)
    z1 = gamma*(r1nne-rcut)
    return phi0*(gamma**3*fzppp(z)/fz(z1))*cutoff(r,rcut)

# define phibar and derivatives w.r.t. r1nn at equilibrium
def phibar(r1nn,phi0,beta,r1nne,rcut):
    rs = zetas*r1nn
    phis = phi(rs,phi0,beta,r1nne,rcut) 
    return np.dot(Zs,phis)

def phibare(phi0,beta,r1nne,rcut):
    rs = zetas*r1nne
    phis = phi(rs,phi0,beta,r1nne,rcut) 
    return np.dot(Zs,phis)

def phibarpe(phi0,beta,r1nne,rcut):
    rs = zetas*r1nne
    phips = phip(rs,phi0,beta,r1nne,rcut) 
    return np.dot(Zs,phips*zetas)

def phibarppe(phi0,beta,r1nne,rcut):
    rs = zetas*r1nne
    phipps = phipp(rs,phi0,beta,r1nne,rcut) 
    return np.dot(Zs,phipps*zetas**2)

def phibarpppe(phi0,beta,r1nne,rcut):
    rs = zetas*r1nne
    phippps = phippp(rs,phi0,beta,r1nne,rcut) 
    return np.dot(Zs,phippps*zetas**3)

# define Ue=U(r1nne), Upe=U', Uppe=U'', Uppe=U'''
def Ue(r1nne,Ece,Be):
    return -Ece

def Uppe(r1nne,Ece,Be):
    return 9*Be*r1nne/math.sqrt(2.)

def Upppe(r1nne,Ece,Be):
    return -27*math.sqrt( math.sqrt(2.)*Be**3 * r1nne**3/Ece )

# coefficients in embedding function
def F0(r1nne,Ece,Be,phi0,rho0,beta,rcut):
    return Ue(r1nne,Ece,Be) - phibare(phi0,beta,r1nne,rcut)/2.

def F1(r1nne,Ece,Be,phi0,rho0,beta,rcut):
    return -phibarpe(phi0,beta,r1nne,rcut)/(2.*rhobarpe(rho0,beta,r1nne,rcut))

def F2(r1nne,Ece,Be,phi0,rho0,beta,rcut):
    Fpe = F1(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    U2 = Uppe(r1nne,Ece,Be)
    phi2 = phibarppe(phi0,beta,r1nne,rcut)
    rho1 = rhobarpe(rho0,beta,r1nne,rcut)
    rho2 = rhobarppe(rho0,beta,r1nne,rcut)
    return (U2-phi2/2.-Fpe*rho2)/rho1**2

def F3(r1nne,Ece,Be,phi0,rho0,beta,rcut):
    Fpe = F1(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fppe = F2(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    U3 = Upppe(r1nne,Ece,Be)
    phi3 = phibarpppe(phi0,beta,r1nne,rcut)
    rho1 = rhobarpe(rho0,beta,r1nne,rcut)
    rho2 = rhobarppe(rho0,beta,r1nne,rcut)
    rho3 = rhobarpppe(rho0,beta,r1nne,rcut)
    return (U3 - phi3/2. - Fpe*rho3 - 3.*Fppe*rho1*rho2)/rho1**3
    
def F4(r1nne,Ece,Be,phi0,rho0,beta,rcut):
    Fc0 = F0(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fc1 = F1(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fc2 = F2(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fc3 = F3(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    rhobe = rhobare(rho0,beta,r1nne,rcut)
    return -24.*(Fc0-Fc1*rhobe+Fc2*rhobe**2/2.-Fc3*rhobe**3/6.)/rhobe**4

def F(rhobar,phi0,rho0,beta,r1nne,Be,Ece,rcut):
    rhobe = rhobare(rho0,beta,r1nne,rcut)
    Fd0 = F0(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fd1 = F1(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fd2 = F2(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fd3 = F3(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fd4 = F4(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    drhobar = rhobar - rhobe
    return Fd0 + Fd1*drhobar + Fd2*drhobar**2/2. + Fd3*drhobar**3/6. + Fd4*drhobar**4/24.

# define criteria for bounds on parameters
# check params against those criteria and flag if they fail

def paramOK(phi0,rho0,beta,r1nne,Be,Ece,rcut):
    Fe0 = F0(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fe1 = F1(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fe2 = F2(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fe3 = F3(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    Fe4 = F4(r1nne,Ece,Be,phi0,rho0,beta,rcut)
    
    crit1 = np.heaviside( Fe2, 0.5 )    # F2>0
    crit2 = np.heaviside( Fe4, 0.5 )    # F4>0
    crit3a = np.heaviside( 2*Fe2*Fe4-Fe3**2, 0.5 ) # 2*F2*F4-F3^2>0
    crit3b = np.heaviside( Fe3, 0.5)  #  F3>0
    crit3 = 0
    if crit3a==1 or crit3b==1:    #  3 = Or(3A,3B)
        crit3 = 1

    rcutmax = math.sqrt(nshellmax+1)*r1nne  # max is figured from settings at top for FCC only 
    crit4 = np.heaviside( rcutmax - rcut, 0.5 )   

    crits = crit1*crit2*crit3*crit4

    if crits != 1:
        print(" parameters failed criteria ")
        print(" crits = ",crit1,crit2,crit3,crit4)
        print(" F0 = ",F0)
        print(" F1 = ",F1)
        print(" F2 = ",F2)
        print(" F3 = ",F3)
        print(" F4 = ",F4)
        print(" rcutmax = ",rcutmax)
        print(" rcut = ",rcut)
        

    return crits

