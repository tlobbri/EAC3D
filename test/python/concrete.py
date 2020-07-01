# -*- coding: utf-8 -*-
"""
Created on Sun Mar 04 17:05:16 2018

@author: tlo
"""

import numpy as np
import math

from scipy.interpolate import splrep, splev

class concrete() :
    def __init__(self):
        self.Cfck = None
        self.cimentType = None
        self.model_type = None
        
        self.fcm_28 = None
        self.fctm_28 = None
        self.Ecm_28 = None
        self.h0 = None
        self.RH = None
        
    
        
    def set_cimentType(self, cimentType='N'):
        self.cimentType = cimentType
        
    def set_modelCode(self, md='EN1992'):
        self.model_type = md
        
    def set_RH(self, RH=60.):
        self.RH = RH
    def set_h0(self, h0):
        self.h0 = h0
        
        
class concrete_MC2010(concrete):
    def __init__(self):
        concrete.__init__(self)
        
        
    def __alpha_ds1(self):
        if self.cimentType == 'S':
            return 3.
        elif self.cimentType == 'N':
            return 4.
        elif self.cimentType == 'R':
            return 6.
        else : 
            return None
        
    def alpha_ds2(self):
        if self.cimentType == 'S':
            return 0.013
        elif self.cimentType == 'N':
            return 0.012
        elif self.cimentType == 'R':
            return 0.012
        else : 
            return None
        
    def set_Cfck(self, Cfck):
        self.Cfck = Cfck
        self.fcm_28 = self.Cfck +8
        if self.Cfck < 50.:
            self.fctm_28 = 0.3*math.pow(self.Cfck, 2./3.)
        else :
            self.fctm_28 = 2.12 * np.log(1.+(self.fcm_28/10.))
            
        alpha_E = 1.0#quartzite
        
        self.Eci_28 = 21.5*math.pow(self.fcm_28/10., 1./3.)*alpha_E
        
        coef = 0.8 + 0.2*self.fcm_28/88.
        if coef > 1.0 : 
            coef = 1.
        self.Ec_28 =  coef* self.Eci_28
        
        
    def fcm(self, time=28.):
        return self.Beta_cc(time)*self.fcm_28
    
    def Gf(self, time=28.) : 
        return 73*self.fcm(time)**0.18
    def __s_Bcc(self):
        
        if self.Cfck > 60.:
            return 0.20
        
        if self.cimentType == 'S':
            return 0.38
        elif self.cimentType == 'N':
            return 0.25
        elif self.cimentType == 'R':
            return 0.20
        else : 
            return None
        
        
    def Beta_cc(self, time):
        return np.exp(self.__s_Bcc() * (1.-(28./time)**0.5))
    
    
    def Eci(self, time):
        Beta_E = self.Beta_cc(time)**0.5
        
        return Beta_E * self.Eci_28
    
    
   
    
    def __Beta_ds(self, time, ts):
        coef = np.ones(time.shape)
        myID  = np.where(time <= ts)
        coef[myID[0]] = 0.0
        
        return coef*((time-ts)/((time-ts)+0.035*self.h0**2))**0.5
        
    
    
    
    def dryingShrinkage(self, time, ts):
        
        Beta_s1 = (35./self.fcm_28)**0.1
        if Beta_s1 > 1.0:
            Beta_s1 = 1.
            
        if (self.RH >= 40) and (self.RH <99.*Beta_s1):
            Beta_RH = 1.55*(1.-math.pow(self.RH/100., 3.))
        elif (self.RH >= 99.*Beta_s1):
            Beta_RH = -0.25
        else:
            Beta_RH = None
        
        
        eps_cds0 = ((220.+110.*self.__alpha_ds1())*np.exp(-self.alpha_ds2()*self.fcm_28))
        
        return eps_cds0*Beta_RH * self.__Beta_ds(time, ts) 
    
    
    def Beta_bs(self, time):
        return 1.-np.exp(-0.2*time**0.5)
    
    def alpha_bs(self):
        if self.cimentType == 'S':
            return 800.
        elif self.cimentType == 'N':
            return 700.
        elif self.cimentType == 'R':
            return 600.
        else : 
            return None
    def autoShrinkage(self, time):
        
        
        eps_cbs0 = self.alpha_bs()*(0.1*self.fcm_28/(6.+0.1*self.fcm_28))**2.5 
        return self.Beta_bs(time)*eps_cbs0
    
    #NOT FOUND IN MC2010!!! take the same than in EC1992
    def fctm(self, time):
        alpha = np.ones(time.shape)
        myID  = np.where(time >= 28)
        alpha[myID[0]] = 2./3.
        return (self.Beta_cc(time)**alpha)*self.fctm_28
    
    def totalShrinkage(self, time, ts):
        return self.autoShrinkage(time) + self.dryingShrinkage(time, ts)
    
    
class concrete_EN1992(concrete):
    def __init__(self):
        concrete.__init__(self)
        
        x = np.array([0, 100., 200., 300., 500., 600.])
        y = np.array([1., 1., 0.85, 0.75, 0.7, 0.7])
        
        self.tick_kh = splrep(x,y,k=1)
        
        
    def set_Cfck(self, Cfck):
        self.Cfck = Cfck
        self.fcm_28 = self.Cfck +8
        if self.Cfck < 50.:
            self.fctm_28 = 0.3*math.pow(self.Cfck, 2./3.)
        else :
            self.fctm_28 = 2.12 * np.log(1.+(self.fcm_28/10.))
        self.Ecm_28 = 22.*math.pow(self.fcm_28/10., 0.3)*1000.
        
        
    def __alpha_ds1(self):
        if self.cimentType == 'S':
            return 3.
        elif self.cimentType == 'N':
            return 4.
        elif self.cimentType == 'R':
            return 6.
        else : 
            return None
        
    def __alpha_ds2(self):
        if self.cimentType == 'S':
            return 0.13
        elif self.cimentType == 'N':
            return 0.12
        elif self.cimentType == 'R':
            return 0.11
        else : 
            return None
    def __s_Bcc(self):
        if self.cimentType == 'S':
            return 0.38
        elif self.cimentType == 'N':
            return 0.25
        elif self.cimentType == 'R':
            return 0.20
        else : 
            return None
        
    def __Beta_cc(self, time):
        return np.exp(self.__s_Bcc() * (1.-(28./np.array(time))**0.5))
    
    def kh(self):
        return splev(self.h0, self.tick_kh)
    
    
    def __Beta_ds(self, time, ts):
        coef = np.ones(time.shape)
        myID  = np.where(time < ts)
        coef[myID[0]] = 0.0
        
        return ((time-ts)/((time-ts)+0.04*math.sqrt(self.h0**3)))*coef
        
    
    def dryingShrinkage(self, time, ts):
            
        Beta_RH = 1.55*(1.-math.pow(self.RH/100., 3.))
        eps_cd0 = 0.85*((220.+110.*self.__alpha_ds1())*np.exp(-self.__alpha_ds2()*self.fcm_28/10.))*Beta_RH
        
        return eps_cd0*self.__Beta_ds(time, ts)*self.kh()
    
    
    def __Beta_as(self, time):
        return 1.-np.exp(-0.2*(np.array(time))**0.5)
    
    def autoShrinkage(self, time):
        eps_ca = 2.5*(self.Cfck-10)
        return self.__Beta_as(time)*eps_ca
    
    def fcm(self, time):
        return self.__Beta_cc(time)*self.fcm_28
    
    def Gf(self, time) : 
        return 73*self.fcm(time)**0.18  #[N/m]
    
    def fctm(self, time):
        time = np.array(time)
        alpha = np.ones(time.shape)
        myID  = np.where(time >= 28)
        alpha[myID[0]] = 2./3.
        return (self.__Beta_cc(time)**alpha)*self.fctm_28
    
    def Ecm(self, time):
        return (self.fcm(time)/self.fcm_28)**0.3 *self.Ecm_28


    def totalShrinkage(self, time, ts):
        return self.autoShrinkage(time) + self.dryingShrinkage(time, ts)
    
    def __BetaCreep_c(self,time, t0):
        
        if self.fcm_28 <=35.:
            BetaH = 1.5*(1+(0.012*self.RH)**18.)*self.h0 + 250.
            if BetaH>1500.:
                BetaH = 1500.
        else:
            
            alpha3 = (35./self.fcm_28)**0.5
            BetaH = 1.5*(1+(0.012*self.RH)**18.)*self.h0 + 250.*alpha3
            
            if BetaH>(1500.*alpha3):
                BetaH = 1500.*alpha3
        
        return ((time-t0)/(BetaH + time-t0))**0.3
    def __phiCreep0(self, t0):
        
        if self.fcm <= 35.:
            phiCreepRH = 1. + (1-self.RH/100.)/(0.1*self.h0**(1./3.))
        else:
            alpha1 = (35./self.fcm_28)**0.7
            alpha2 = (35./self.fcm_28)** 0.2
            phiCreepRH = (1. + (1-self.RH/100.)/(0.1*self.h0**( 1./3.))*alpha1)*alpha2
        BetaCreep_fcm = 16.8/math.sqrt(self.fcm_28)
        BetaCreep_t0 = 1./(0.1+t0**0.2)
        
        return phiCreepRH*BetaCreep_fcm*BetaCreep_t0
    
    
    def totalphiCreep(self, time, t0):
        return self.__phiCreep0(t0)*self.__BetaCreep_c(time, t0)