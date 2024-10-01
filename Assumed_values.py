# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 15:39:54 2021

@author: Moges Retta
List of assumed values that are constant for species

ei              # facilitation factor for membrane or stroma,
                  unitless 
fi              # fraction of the diffusive path length of component and its 
                  thickness, 
                  m/m
Gmem            # Permeability of the cell wall
                  m/s

Genv            # Permeability of the chloroplast envelope, 
                  m/s
KmC             # Michaelis–Menten constant of Rubisco for CO2, 
                  µbar CO2

KmO             # Michaelis–Menten constant of Rubisco for O, 
                  mbar
peff_i          # Effective porosity of component i, unitless

ti              # Weighted average thickness of a mesophyll component i in the 
                  palisade and the spongy parenchyma, m
                  
Sco             # Relative CO2/O2 specificity factor of Rubisco, 
                              mbar O2 µbar−1 CO2                      

gamma_tissue    # Curvature factor of a certain tissue 
                  (either palisade parenchyma or spongy parenchyma)
                  
zeta_i          # Reduction factor of the diffusion coefficient of CO2 
                  relative to DCO2,water in component i due to the higher 
                  viscosity of i
                  
"""
class Assumed_values:
    KmC = 267             
    KmO = 164  
    sco = 2.86;

    # cell wall
    f_wall = 1.0;
    peff_wall = 0.3;
    zeta_wall = 1.0;
    f_wall_spo = 1.0;
    peff_wall_spo = 0.3;
    zeta_wall_spo = 1.0;
    e_cw = 0; 

    
    # pm
    Gmem = 3.5*10**-3; # facilitation, x5
    zeta_pm = 1.0;
    peff_pm= 1;
    f_pm = 1.0;
    t_pm = 1;
    e_pm = 0; 
    
    # chl. env
    f_env= 1.0;
    t_env = 1;    
    peff_env = 1;
    zeta_env = 1.0;
    Genv = Gmem/2
    e_ce = e_pm;

   
    # cytosol
    f_cytsol = 1.0; # facilitation, /5
    zeta_cytsol = 0.50;
    peff_cytsol = 1;
    f_cytsol_spo = 1.0; # facilitation, /5
    zeta_cytsol_spo = 0.50;
    peff_cytsol_spo = 1;
    e_cy = 0; 

    # stroma  
    f_stroma = 0.25; # facilitation, /10        
    peff_stroma = 1;
    zeta_stroma = 0.5;
    f_stroma_spo = 0.5; # facilitation, /10        
    peff_stroma_spo = 1;
    zeta_stroma_spo = 0.5;
    e_st = 0; 

    
    # curvature factors
    gamma_tissue_pali = 1             
    gamma_tissue_spongy= 1 # 3D values were used
    
    # intercellular airspace
    turtosity = 1.57 # recent 3-D value 1.16, expression (porosity)**-0.18=1.28
    factor =4 # amphistomatous leaves has effective length of 1/4th Lmesophyll
    path_lengthening = 1; # 1.42 for C3 and 1.78 for CAM
    connectivity = 1; # 0.94 for c3 and 0.88 for C4
    
    def __init__(self):
        self.values = []
      
        
    def get_kinetic_constants(self):
        return [self.KmC,self.KmO,self.sco]
    
    
    def get_curvature(self):
        return[self.gamma_tissue_pali,self.gamma_tissue_spongy]
    
    
    def get_stroma(self):
        return [self.f_stroma,self.peff_stroma,self.zeta_stroma,self.e_st]
    
    
    def get_cytosol(self):
        return [self.f_cytsol,self.peff_cytsol,self.zeta_cytsol,self.e_cy]
    
    
    def get_cell_wall(self):
        return [self.f_wall,self.peff_wall,self.zeta_wall,self.e_cw]

    def get_stroma_spo(self):
        return [self.f_stroma_spo,self.peff_stroma_spo,self.zeta_stroma_spo,self.e_st]
    
    
    def get_cytosol_spo(self):
        return [self.f_cytsol_spo,self.peff_cytsol_spo,self.zeta_cytsol_spo,self.e_cy]
    
    
    def get_cell_wall_spo(self):
        return [self.f_wall_spo,self.peff_wall_spo,self.zeta_wall_spo,self.e_cw]


    def get_plasmalemma(self):
        return [self.f_pm,self.peff_pm,self.t_pm,self.zeta_pm,self.Gmem,self.e_pm]


    def get_chloroplast_envelope(self):
        return [self.f_env,self.peff_env,self.t_env,self.zeta_env,self.Genv,self.e_ce]
    
    def get_turtosity_ias(self):
        return self.turtosity