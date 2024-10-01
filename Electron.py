# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 15:37:49 2021

@author: Moges Retta
Jmax            # Maximum rate of electron transport through Photosystem II 
                  at saturating light, 
                  µmol e− m−2 leaf s−1
                  
teta            # Convexity factor of the response of J to Iinc                  

                  
alpha2LL        # Quantum yield of electron transport through Photosystem II
                  under strictly electron-transport-limiting conditions on the
                  basis of light absorbed by both Photosystem I and Photosystem
                  II,
                  mol e− mol−1 photon

phi2            # Quantum yield of electron transport through Photosystem II, 
                  mol e− mol−1 photon
                  
k2LL            # Conversion factor of incident irradiance into electron 
                  transport under electron-transport-limited conditions, 
                  mol e− mol−1 photon

s               # Slope of the assumed linear relationship between AN and 1
                # 1/4*Iinc*Phi2 under strictly electron-transport-limited conditions
                                    
"""  
import math

class Electron:
    def __init__(self,Jmax,k2LL,teta,phi2,Iinc):
        self.set_Jmax(Jmax)
        self.set_k2LL(k2LL)
        self.set_teta(teta)
        self.set_phi2(phi2)
        self.set_Iinc(Iinc)
        
    def set_Jmax(self,Jmax):
        if Jmax>=0:
            self.Jmax = Jmax
    def set_k2LL(self,k2LL):
        if k2LL>=0:
            self.k2LL = k2LL     
    def set_teta(self,teta):
        if teta>=0:
            self.teta = teta
    def set_phi2(self,phi2):
        if phi2>=0:
            self.phi2 = phi2
    def set_Iinc(self,Iinc):
        if Iinc>=0:
            self.Iinc = Iinc
    
    def get_Jmax(self):
        return self.Jmax
    def get_k2LL(self):
        return self.k2LL
    def get_teta(self):
        return self.teta
    def get_phi2(self):
        return self.phi2
    def get_Iinc(self):
        return self.Iinc
    
    def calculated_J(self):
        A = self.get_k2LL()*self.get_Iinc()+self.get_Jmax()
        B = (4*self.get_teta()*self.get_Jmax()*self.get_k2LL()*self.get_Iinc())
        return (A-math.sqrt(A**2-B))/(2*self.get_teta())