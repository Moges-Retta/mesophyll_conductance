# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 16:00:07 2021

@author: Moges
"""
import pandas as pd
from Assumed_values import Assumed_values
import numpy as np

class Data:
    
    """ Constants and measurement data for anatomy and photosynthesis kinetics
        
         ----------------------------------------------------------------------
        DCO2_i          # Diffusion coefficient of CO2 in component i, 
                          m2/s
        
        DCO2,water      # Diffusion coefficient of CO2 in water, 
                          m2/s
        e_i              # facilitation factor, unitless     
        
        fi              # fraction of the diffusive path length of component and its 
                          thickness, 
                          m/m
        
        fpal            # fraction of the exposed mesophyll surface area that belongs 
                          to the palisade
                          m/m
        
        F               # Rate of photorespiratory CO2 release, 
                          µmol CO2 m−2 leaf s−1
        
        Gmem            # Permeability of the cell wall
                          m/s
        
        Genv            # Permeability of the chloroplast envelope, 
                          m/s
        
        H               # Henry’s law constant for CO2,
                          Pa m−3 mol−1
        
        Iinc            # Irradiance incident at the leaf surface, 
                          µmol/m2/s
        
        J               # Rate of electron transport through Photosystem II, 
                          µmole− m−2 leaf s−1
        
        Jmax            # Maximum rate of electron transport through Photosystem II 
                          at saturating light, 
                          µmol e− m−2 leaf s−1
                          
        KmC             # Michaelis–Menten constant of Rubisco for CO2, 
                          µbar CO2
        
        KmO             # Michaelis–Menten constant of Rubisco for O, 
                          mbar
        
        Li              # Diffusion path length of component i, 
                          µm
        
        Lm_L_tissue     # Fraction of exposed mesophyll length relative to the width 
                          of the section at one side of the leaf in a certain tissue 
                          (either palisade parenchyma or spongy, 
                          m/m
                          
        O               # O2 partial pressure, 
                          mbar
        
        peff_i          # Effective porosity of component i, unitless
        
        q               # Power in the power law that describes the empirical 
                          relationship between Ci and Iinc
                          
        ri              # Resistance for CO2 transport of component i in the 
                          mesophyll, 
                          m2 leaf s bar CO2 mol−1 CO2
                          
        rchi            # Lumped resistance for CO2 transport of the chloroplast 
                          envelope and the stroma, and half the resistance of the 
                          cytosol, 
                          m2 leaf s bar CO2 mol−1 CO2
                          
        rdiff           # Total resistance for CO2 transport of the physical barriers 
                          in the mesophyll, 
                          m2 leaf s bar CO2 mol−1 CO2
                          
        rm              # Apparent mesophyll resistance, 
                          m2 leaf s bar CO2 mol−1 CO2
        
        rwp             # Lumped resistance for CO2 transport of the cell wall, 
                          the plasmamembrane, and half the resistance of the cytosol,
                          m2 leaf s bar CO2 mol−1 CO2
                          
        R               # Universal gas constant, 
                          JK−1 mol−1
                          
        Rd              # Rate of mitochondrial respiration in the light, 
                          µmol CO2 m−2 leaf s−1
                          
        Ri              # Resistance for CO2 transport of component in the mesophyll,
                          m2 leaf s bar CO2 mol−1 CO2
                          
        s               # Slope of the assumed linear relationship between AN and 1
                        # 1/4*Iinc*Phi2 under strictly electron-transport-limited conditions
                        
        Sc_S            # Fraction of the exposed chloroplast surface area of the 
                          palisade parenchyma and the spongy parenchyma relative to 
                          leaf surface area at one one side of the leaf, m/m
                          
        Sc_Sm           # Fraction of the exposed chloroplast surface area of the palisade
                        #   parenchyma and the spongy parenchyma relative to the 
                        exposed mesophyll surface area of these tissues
                        
        Sco             # Relative CO2/O2 specificity factor of Rubisco, 
                          mbar O2 µbar−1 CO2
                          
        Sm_S            # Fraction of the exposed mesophyll surface area of the
                          palisade parenchyma and the spongy parenchyma relative to 
                          leaf surface area at one side of the leaf, m/m
                          
        ti              # Weighted average thickness of a mesophyll component i in the 
                          palisade and the spongy parenchyma, m
                          
        T               # Temperature, °K
        
        Tp              # Rate of triose phosphate utilization, 
                          µmol phosphate m−2 leaf s−1
                          
        alpha2LL        # Quantum yield of electron transport through Photosystem II
                          under strictly electron-transport-limiting conditions on the
                          basis of light absorbed by both Photosystem I and Photosystem
                          II,
                          mol e− mol−1 photon
                          
        gamma_tissue    # Curvature factor of a certain tissue 
                          (either palisade parenchyma or spongy parenchyma)
                          
        gamma_star      # CO2 compensation point, 
                          µbar CO2
        teta            # Convexity factor of the response of J to Iinc
        
        zeta_i          # Reduction factor of the diffusion coefficient of CO2 
                          relative to DCO2,water in component i due to the higher 
                          viscosity of i
                          
        k2LL            # Conversion factor of incident irradiance into electron 
                          transport under electron-transport-limited conditions, 
                          mol e− mol−1 photon
                          
        phi2            # Quantum yield of electron transport through Photosystem II, 
                          mol e− mol−1 photon
                          
        ω               # Ratio of rchl to rdiff
        -------------------------------------------------------------------------------
        """
    data = pd.read_excel ('Parameters.xlsx',sheet_name='Leaf_anatomy') 
    FORMAT = ['Plant','Treatment', 'Replicate','Tissue','Component', 'Thickness', 'Lm L','Lc Lm','Mesophyll_thickness','Porosity']
    leaf_anatomy = data[FORMAT]
    data = pd.read_excel ('Parameters.xlsx',sheet_name='Kinetic_parameters') 
    FORMAT = ['Plant','Treatment','Replicate','Vcmax','Rd', 'Sco', 'Tp','Jmax','k2LL','theta','gm','alphaS','sigma']
    kinetic_constants = data[FORMAT]
    assumed_values = Assumed_values()
    [gamma_tissue_pali,gamma_tissue_spongy]= assumed_values.get_curvature()
    [KmC,KmO,sco]= assumed_values.get_kinetic_constants()
    replicates = kinetic_constants['Replicate'].unique()
    
    def __init__(self,plant):
            self.plant = plant
          
             
    def get_data_frame(self,DF):
        plant_name = self.plant.get_name()
        treatment = self.plant.get_treatment()
        df = DF[DF['Plant']==plant_name]
        df = df[df['Treatment']==treatment]
        return df
    
    # Leaf anatomical properties measured from light microscopy
    # Curvature factors, gamma_tissue are assumed
    def get_leaf_anatomy(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Tissue','Replicate','Lm L','Lc Lm']]
        df = df[df['Replicate']==replicate]

        pali = df[df['Tissue']=='Palisade']
        spon = df[df['Tissue']=='Spongy']

        Lm_L_tissue_pali = pali[['Lm L']].values
        Lm_L_tissue_spongy = spon[['Lm L']].values
        Lc_Lm_tissue_pali = pali[['Lc Lm']].values
        Lc_Lm_tissue_spongy = spon[['Lc Lm']].values

        df.loc[:,'gamma_tissue_pali']=self.gamma_tissue_pali
        df.loc[:,'gamma_tissue_spongy']=self.gamma_tissue_spongy
        
        f_pal = (Lc_Lm_tissue_pali*Lm_L_tissue_pali/(Lc_Lm_tissue_pali*Lm_L_tissue_pali+Lc_Lm_tissue_spongy*Lm_L_tissue_spongy));
        df.loc[:,'f_pal']=f_pal[0][0]

        return df
    

    def get_mesophyll_thickness(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Replicate','Mesophyll_thickness']]
        df = df[df['Replicate']==replicate]

        return df['Mesophyll_thickness'].values[0]
    
    
    def get_leaf_porosity(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Replicate','Porosity']]
        df = df[df['Replicate']==replicate]
        return df['Porosity'].values[0]
    
    # FvCB model kinetic constants for Rubisco        
    def get_rubisco_kinetics(self):
        df = self.get_data_frame(self.kinetic_constants)
        df = df[['Replicate','Vcmax','Rd','Sco','Tp','gm','alphaS']]
        df.loc[:,'KmC']=self.KmC
        df.loc[:,'KmO']=self.KmO
        return df
    
    
    # Parameters of non-rectangular hyperbolic equation for J        
    def get_electron(self):
        df = self.get_data_frame(self.kinetic_constants)
        df = df[['Replicate','Jmax','k2LL','theta']]
        return df
    
    # Get sigma of gm model, Yin        
    def get_sigma(self):
        df = self.get_data_frame(self.kinetic_constants)
        return df.loc[:,'sigma'].unique()
    
    # Cell wall of spongy cells
    def get_spongy_wall(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Tissue','Replicate','Component','Thickness']]
        spon = df[df['Tissue']=='Spongy']
        spon = spon[spon['Component']=='Cell wall']
        spon = spon[spon['Replicate']==replicate]

        ti_spo =  spon['Thickness'].values
        [fi,peff_i,zeta_i,e_i] = self.assumed_values.get_cell_wall_spo()

        return [fi,ti_spo,peff_i,zeta_i,e_i]
    
    
    # Cytsol of spongy cells
    def get_spongy_cytosol(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Tissue','Replicate','Component','Thickness']]
        spon = df[df['Tissue']=='Spongy']
        spon = spon[spon['Component']=='Cytosol']
        spon = spon[spon['Replicate']==replicate]

        ti_spo =  spon['Thickness'].values
        [fi,peff_i,zeta_i,e_i] = self.assumed_values.get_cytosol_spo()

        return [fi,ti_spo,peff_i,zeta_i,e_i]
    
    
    # Stroma of spongy cells
    def get_spongy_stroma(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Tissue','Replicate','Component','Thickness']]
        spon = df[df['Tissue']=='Spongy']
        spon = spon[spon['Component']=='Stroma']
        spon = spon[spon['Replicate']==replicate]

        ti_spo =  spon['Thickness'].values[0]
        [fi,peff_i,zeta_i,e_i] = self.assumed_values.get_stroma_spo()
        
        return [fi,ti_spo,peff_i,zeta_i,e_i]
    
    
    # Plasma membrane , resistance modeled as 1/Gmem
    def get_plasmallema(self):
        [fi,peff_i,ti,zeta_i,Gmem,e_i] = self.assumed_values.get_plasmalemma()          
        return [fi,ti,peff_i,zeta_i,Gmem,e_i]
    
    # Chloroplast envelop , resistance modeled as 1/Genv
    def get_chlenvelop(self):
        [fi,peff_i,ti,zeta_i,Genv,e_i] = self.assumed_values.get_chloroplast_envelope()          
        return [fi,ti,peff_i,zeta_i,Genv,e_i]
    
    
    # Cell wall of palisade cells    
    def get_palisade_wall(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Tissue','Replicate','Component','Thickness']]

        pali = df[df['Tissue']=='Palisade']
        pali = pali[pali['Component']=='Cell wall']
        pali = pali[pali['Replicate']==replicate]

        ti_pal =  pali['Thickness'].values[0]
        [fi,peff_i,zeta_i,e_i] = self.assumed_values.get_cell_wall()
        
        return [fi,ti_pal,peff_i,zeta_i,e_i]
    
    
    # Cytsol of palisade cells
    def get_palisade_cytosol(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Tissue','Replicate','Component','Thickness']]
        pali = df[df['Tissue']=='Palisade']
        pali = pali[pali['Component']=='Cytosol']
        pali = pali[pali['Replicate']==replicate]
        
        ti_pal =  pali['Thickness'].values[0]
        [fi,peff_i,zeta_i,e_i] = self.assumed_values.get_cytosol()
 
        return [fi,ti_pal,peff_i,zeta_i,e_i]
    
    
    # Stroma of spongy cells
    def get_palisade_stroma(self,replicate):
        df = self.get_data_frame(self.leaf_anatomy)
        df = df[['Tissue','Replicate','Component','Thickness']]
        pali = df[df['Tissue']=='Palisade']
        pali = pali[pali['Component']=='Stroma']
        pali = pali[pali['Replicate']==replicate]

        ti_pal =  pali['Thickness'].values[0]
        [fi,peff_i,zeta_i,e_i] = self.assumed_values.get_stroma()

        return [fi,ti_pal,peff_i,zeta_i,e_i]