# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 14:05:30 2019

@author: tlo
"""

import sys,os
from EAC3D_module.EAC3D_module import *
import time
import numpy as np
from scipy.interpolate import splev, splrep
import concrete
import math

import matplotlib.pyplot as plt

def alpha_func(RH, cType, fcm, h0):
    
    """def kh(h0):
        x = np.array([0, 100., 200., 300., 500., 600.])
        y = np.array([1., 1., 0.85, 0.75, 0.7, 0.7])
        tick = splrep(x,y,k=1)
        
        return splev(h0, tick)    
    """
    ########################################
    #  IF scipy and numpy not installed    #
    #   use the func below  instead        #
    #    --> h0 = 100 [mm] only!!!!        #
    ########################################
    
    def kh(h0):
        return 1.0
    
    def alpha_ds1(cType):
        if cType == 'S':
            return 3.
        elif cType == 'N':
            return 4.
        elif cType == 'R':
            return 6.
        else : 
            return None
        
    def alpha_ds2(cType):
        if cType == 'S':
            return 0.13
        elif cType == 'N':
            return 0.12
        elif cType == 'R':
            return 0.11
        else : 
            return None

    Beta_RH = -1.55*(1.-math.pow(RH, 3.))
    
    eps_cd0 = 0.85*((220.+110.*alpha_ds1(cType))*math.exp(-alpha_ds2(cType)*fcm/10.))*Beta_RH*1e-6
    drying_shrinkage = eps_cd0*kh(h0)

    
    return drying_shrinkage

h0  = 40.
cType = 'R'
fck = 33

myConcrete = concrete.concrete_EN1992()
myConcrete.set_Cfck(fck)
myConcrete.set_RH(60.)
myConcrete.set_h0(h0)
myConcrete.set_cimentType(cType)
fcm = myConcrete.fcm_28

timeVect = np.linspace(0,200,200)
shrinkage = myConcrete.dryingShrinkage(timeVect, 7)

plt.figure()
plt.semilogy(timeVect, shrinkage)

x = np.linspace(0, 1.1,500)
y = np.zeros(x.shape)
for i in range(len(y)):
    y[i] = alpha_func(x[i], cType, fcm, h0)
tck = splrep(x,y)

plt.figure()
plt.semilogy(x,(y-y[0])*1e6)   

#sys.exit(-1)
def alpha_func_der(RH):
    return splev(RH, tck, der=0)
    
    #return (splev(RH, tck) - splev(0,tck)) / (RH)
   

def bnd_fixed( y,  z,  time, T):
    return 1.0

def bnd_conv( y,  z,  time, H):
    #hD = 1.5e-9#1/s
    #hD = 1e-9*200000.0#300#1e-7*20000000.0#1/s
    hD =1e-09#1e-9#0.001431242#*2e-6
    #if H<60.:
    #    H=60.
    return (0.6-H)*hD*10
    #if H<0. :  H = 0.1
    #return math.log(0.6/H)*hD#*DH(H)


def bnd_isolated(x,y,  time, T):
    return 0.0



def DH(H):
    fcm = 35+8
    D10 = 1.135e-8
    D1 = D10/(fcm-8.)
    alphaD = 0.074
    Hc = 0.722
    n = 10.86
    
    #H = H/100.

    if (H>=1):
        H = 1
    elif H<0:
        H = 0
    
    D = D1*(alphaD+(1.-alphaD)/(1+((1.-H)/(1.-Hc))**n))
    #Dcopper = 401/(390*8960.)
    
    return D
    #return D10*2e6
def DStrainODH(H):
    return alpha_func_der(H)
    


def Restrain(x,y,z):
    return 0.75

def J_relax(t,tprime):
    
    tfactor = 3600.*24.
    creep = True
    #myConcrete = concrete.concrete_EN1992()
    #myConcrete.set_Cfck(fck)
    #myConcrete.set_RH(60.)
    #myConcrete.set_h0(h0)
    #myConcrete.set_cimentType(cType)
#    return 1./myConcrete.totalphiCreep(t,tprime)
    if creep : 
        return (1./myConcrete.Ecm(tprime/tfactor)+myConcrete.totalphiCreep(t/tfactor,tprime/tfactor)/myConcrete.Ecm_28)
    else:
        return 1./myConcrete.Ecm_28
    

"""
def DT(T):
    return 7.5e-7;
"""
def init_func_T( x,  y,  z):
    return 0.0#50.*z;
"""
def init_func_H( x,  y,  z):
    return 1.0;
"""
Lx = 0.01
Ly = 0.01
Lz = 1.0
nx = 5
ny = 5
nz = 100
    
    
mySolver = EAC3D(Lx, Ly, Lz, nx, ny, nz)
mySolver.set_NBEquations(1);
mySolver.saveMesh()


#DEFINE INITIAL TEMPERATURE FIELD AND ITS BND COND.
mySolver.init_field_py(init_func_T, TID)

###############
#   XS
##############
mySolver.set_T_BNDXS_py(bnd_isolated, BND_NEU)
#mySolver.set_T_BNDXS_py(bnd_fixed, BND_DIR)


##################
#   XF
#################
#mySolver.set_T_BNDXF_py(bnd_conv, BND_NEU)
mySolver.set_T_BNDXF_py(bnd_isolated, BND_NEU)
#mySolver.set_T_BNDXF_py(bnd_fixed, BND_DIR)

#####################
#     YS
#####################
mySolver.set_T_BNDYS_py(bnd_isolated, BND_NEU)
#mySolver.set_T_BNDYS_py(bnd_fixed, BND_DIR)

##########################
#       YF
##########################
mySolver.set_T_BNDYF_py(bnd_isolated, BND_NEU)
#mySolver.set_T_BNDYF_py(bnd_conv, BND_NEU)
#mySolver.set_T_BNDYF_py(bnd_fixed, BND_DIR)


########################
#         ZS
######################
mySolver.set_T_BNDZS_py(bnd_isolated, BND_NEU)


################################
#           ZF
################################
#mySolver.set_T_BNDZF_py(bnd_isolated, BND_NEU)
mySolver.set_T_BNDZF_py(bnd_fixed, BND_DIR)
mySolver.set_T_lambda_py(DH)

mySolver.set_DStrainODT_py(DStrainODH)
mySolver.set_Restrain_py(Restrain)
mySolver.set_Relaxation_py(J_relax)

#DEFINE INITIAL H FIELD AND ITS BND COND.
"""
mySolver.init_field_py(init_func_H, HID)
mySolver.set_H_BNDXS_py(bnd_xs, BND_NEU)
mySolver.set_H_BNDXF_py(bnd_xf, BND_DIR)

mySolver.set_H_BNDYS_py(bnd_ys, BND_DIR)
mySolver.set_H_BNDYF_py(bnd_yf, BND_DIR)

mySolver.set_H_BNDZS_py(bnd_zs, BND_DIR)
mySolver.set_H_BNDZF_py(bnd_zf, BND_DIR)
mySolver.set_H_lambda_py(DH)
"""

#DEFINE GLOBAL QUANTITIES

mySolver.set_rho(1.0)
mySolver.set_cp(1.0)

t_vect = [0]
            
nb = 100


for i in range(nb):
    t_vect.append(t_vect[-1]+0.1)

for i in range(nb):
    t_vect.append(t_vect[-1]+1.)
    


for i in range(2*nb):
    t_vect.append(t_vect[-1]+10.)


"""
for i in range(nb):
    t_vect.append(t_vect[-1]+1000)
"""    
t_vect = np.array(t_vect)*3600*24


t1 = time.clock()


for tid in range(1,len(t_vect)):
    t = t_vect[tid]
    dt = mySolver.set_dt(t_vect[tid]-t_vect[tid-1])
    rx = dt*DH(1)/(Lx/(nx-1.))**2
    ry = dt*DH(1)/(Ly/(ny-1.))**2
    rz = dt*DH(1)/(Lz/(nz-1.))**2
    print "tid=%d, dt=%.4f time=%.4f [s] , %.4f [d], r = [%.2f %.2f %.2f] " % (tid, dt, t_vect[tid], t_vect[tid]/(3600.*24), rx, ry, rz)
    
    if ((tid % 10) == 0) :
        datafile = ("data_T_%.4d.dat" % tid)
        mySolver.saveData(datafile, TID);
        print datafile
        
        
        datafile = ("data_Stress_%.4d.dat" % tid)
        mySolver.saveData(datafile, STRESSID);
        print datafile
        
        datafile = ("data_Strain_%.4d.dat" % tid)
        mySolver.saveData(datafile, STRAINID);
        print datafile
        
    if False :
        mySolver.timeMarching(True)
    else:
        mySolver.implicitTimeMarching(1.0, 3)
    #mySolver.computeStrainStress()
t2 = time.clock()

print t2-t1