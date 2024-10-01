# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 15:11:28 2021

@author: Moges Retta


                  
KmC             # Michaelis–Menten constant of Rubisco for CO2, 
                  µbar CO2

KmO             # Michaelis–Menten constant of Rubisco for O, 
                  mbar
Tp              # Rate of triose phosphate utilization, 
                  µmol phosphate m−2 leaf s−1

                
Sco             # Relative CO2/O2 specificity factor of Rubisco, 
                  mbar O2 µbar−1 CO2    
                  
gamma_star      # CO2 compensation point, 
                  µbar CO2
                  
"""  

class Rubisco:
    def __init__(self,Vcmax,KmC,KmO,Rd,Tp,Sco):
        self.set_Vcmax(Vcmax)
        self.set_KmC(KmC)
        self.setKmO(KmO)
        self.set_Rd(Rd)
        self.set_Tp(Tp)
        self.set_Sco(Sco)
        
    def set_Vcmax(self,Vcmax):
        if Vcmax>=0:
            self.Vcmax = Vcmax
            
            
    def set_KmC(self,KmC):
        if KmC>=0:
            self.KmC = KmC   
            
            
    def set_KmO(self,KmO):
        if KmO>=0:
            self.KmO = KmO
            
            
    def set_Rd(self,Rd):
        if Rd>=0:
            self.Rd = Rd
            
            
    def set_Tp(self,Tp):
        if Tp>=0:
            self.Tp = Tp
            
            
    def set_Sco(self,Sco):
        if Sco>=0:
            self.Sco = Sco
            
    
    def get_Vcmax(self):
        return self.Vcmax
    
    
    def get_KmC(self):
        return self.KmC
    
    
    def get_KmO(self):
        return self.KmO
    
    
    def get_Rd(self):
        return self.Rd
    
    
    def get_Tp(self):
        return self.Tp
    
    
    def get_Sco(self):
        return self.Sco