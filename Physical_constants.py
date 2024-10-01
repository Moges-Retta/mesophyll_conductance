# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 12:15:13 2021

@author: Moges Retta
Physical constants
T               # Temperature, °K
                  
R               # Universal gas constant, 
                  JK−1 mol−1
                  
H               # Henry’s law constant for CO2,
                  Pa m^3 mol−1
P               # atmospheric pressure, Pa           
"""
from scipy import constants

class Physical_constants:
        T = constants.zero_Celsius + 25; #°K
        P = constants.atm 
        R = constants.R    
        H = 2941; 
        DCO2_water = 1.79e-9 
        DCO2_gas = 1.51e-5 
            
        def __init__(self):
            self.phyical_constants = [self.T,self.P,self.R,self.H,self.DCO2_water]
    
    
    # Physical constants       
        def get_phyical_constants(self):
            return self.phyical_constants
        
        
        def set_physical_constants(self,T,P,R,H,DCO2_water):
            self.T = T
            self.R = R
            self.P = P
            self.H = H
            self.DCO2_water=DCO2_water
            