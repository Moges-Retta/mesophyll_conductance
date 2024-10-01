# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:11:10 2021
Store leaf anatomical and transport properties
@author: Moges
"""
class Cell_component:
        """
        Created on Mon Jan 18 14:11:10 2021
        Store leaf anatomical and transport properties
        @author: Moges
        Model organelle as having :
        ei              # facilitation factor for diffusion, unitless
        fi              # fraction of the diffusive path length of component and its 
                          thickness, 
                          m/m
        ti              # Weighted average thickness of a mesophyll component i in the 
                          palisade and the spongy parenchyma, m
                          
        peff_i          # Effective porosity of component i, unitless
        
        zeta_i          # Reduction factor of the diffusion coefficient of CO2 
                          relative to DCO2,water in component i due to the higher 
                          viscosity of i
        DCO2,water      # Diffusion coefficient of CO2 in water, 
                          m2/s
        """

        def __init__(self,name,f,t,peff,zeta,DCO2_water,ei):
            self.name=name
            self.set_f(f)
            self.set_t(t)
            self.set_peff(peff)
            self.set_zeta(zeta)
            self.set_DCO2_water(DCO2_water)
            self.set_ei(ei)

            
               
        def set_name(self,name):
            if isinstance(name, str):
                self.name=name
            else:
                print("no st")
                
                
        def set_f(self,f):
            if f>=0:
                self.f=f
                
                
        def set_t(self,t):
            if t>0:
                self.t=t
                
                
        def set_peff(self,peff):
            if peff>0:
                self.peff=peff
                
                
        def set_zeta(self,zeta):
            if zeta>=0:
                self.zeta=zeta
                
                
        def set_DCO2_water(self,DCO2_water):
            if DCO2_water>=0:
                self.DCO2_water=DCO2_water
        
        def set_ei(self,ei):
            if ei>=0:
                self.ei=ei
                
            
        def get_name(self):
            return self.name
        
        
        def get_f(self):
            return self.f
        
        
        def get_t(self):
            return self.t
        
        
        def get_peff(self):
            return self.peff
        
        
        def get_zeta(self):
            return self.zeta
        
        
        def get_DCO2_water(self):
            return self.DCO2_water
        
        def get_ei(self):
            return self.ei
        

        def calculate_component_resistance(self):  
            """ calculate resistance of a cell component """
            return self.get_f()*self.get_t()/(self.get_peff()*
                           self.get_zeta()*self.get_DCO2_water()*(1+self.get_ei()));