# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:06:31 2021                  
#@author: Moges
"""

class Cell:
    """
    Created on Mon Jan 18 14:06:31 2021
    Model cell as having surface and total resistance
    Calculate total resistance of components in a cell type, exposed length of 
    mesophyll and chloroplast per leaf surface
    
    fpal            # raction of the exposed mesophyll surface area that belongs 
                      to the palisade
                      m/m
    Lm_L     # Fraction of exposed mesophyll length relative to the width 
                      of the section at one side of the leaf in a certain tissue 
                      (either palisade parenchyma or spongy, 
                      m/m
    Lc_Lm     # Fraction of exposed chloroplast length relative to exposed mesophyll
                surface  (either palisade parenchyma or spongy, 
                      m/m                  
    gamma_tissue    # Curvature factor of a certain tissue 
                      (either palisade parenchyma or spongy parenchyma)
                      
    #@author: Moges
    """

    def __init__(self,name,Lm_L,Lc_Lm,gamma_tissue,f_pal):
            self.name=name
            self.set_Lm_L(Lm_L)
            self.set_Lc_Lm(Lc_Lm)
            self.set_gamma_tissue(gamma_tissue)
            self.set_f_pal(f_pal)
            self.Cell_components=[]  
            
            
    def set_name(self,name):
        if isinstance(name, str):
            self.name=name


    def set_Lm_L(self,Lm_L):
        if Lm_L>=0:
            self.Lm_L=Lm_L
           
            
    def set_Lc_Lm(self,Lc_Lm):
        if Lc_Lm>=0:
            self.Lc_Lm=Lc_Lm 
        
        
    def set_gamma_tissue(self,gamma_tissue):
        if gamma_tissue>=0:
            self.gamma_tissue=gamma_tissue  
         
            
    def set_f_pal(self,f_pal):
        if f_pal>=0:
            self.f_pal=f_pal   
         
            
    def set_curvature(self,gamma_tissue):
        if gamma_tissue>=0:
            self.gamma_tissue=gamma_tissue  
        
        
    def get_name(self):
        return self.name
    
    
    def get_Lm_L(self):
            return self.Lm_L
        
    def get_Lc_Lm(self):
            return self.Lc_Lm
        
        
    def get_curvature(self):
            return self.gamma_tissue
        
    def get_f_pal(self):
            return self.f_pal 
       
        
    def get_cell_components(self):
        return self.Cell_components
    
    
    def add_organelles(self,cell_component):
        self.Cell_components.append(cell_component)
       
    def remove_organelles(self,cell_component):
        for component in self.Cell_components:
            if component.get_name()==cell_component.get_name():
                self.Cell_components.remove(component)
        
    def calculate_Sm_S(self):
        return self.get_Lm_L()*self.get_curvature()
    
    
    def calculate_Sc_S(self):
        return self.get_Lm_L()*self.get_Lc_Lm()*self.get_curvature()