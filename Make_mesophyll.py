# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 21:25:28 2021

@author: Moges
Make mesophyll tissue from inputs of anatomical and transport properties
The class allows construction of cell component from the inputs data or from
an object given as an input and the type of tissue the object belongs
inputs : inputs.py
cell component :  cell component object
types : string input, spongy or palisade tissue
"""
from Cell import Cell
from Cell_component import Cell_component
from Names import Names
from Physical_constants import Physical_constants

class Make_mesophyll:
    
    def __init__(self,inputs,cell_component,types):
        self.inputs=inputs
        self.cell_component=cell_component
        self.types=types
      
        
    def get_inputs(self):
        return self.inputs

    
    def add_palisade_organells(self,palisade):
        physical_constants = Physical_constants()
        [T,P,R,H,DCO2_water] = physical_constants.get_phyical_constants()
        
        """ Add cell components of palisade cells"""
        [fi,ti_pal,peff_i,zeta_i] = self.inputs.get_palisade_wall()
        cell_wall = Cell_component(Names.CELLWALL,fi,ti_pal,peff_i,zeta_i,DCO2_water);
        palisade.add_organelles(cell_wall)
        
        # plasma membrane
        [fi,ti,peff_i,zeta_i,Gmem,e_i]=self.inputs.get_plasmallema()
        plasma_membrane = Cell_component(Names.PLASMALLEMA,fi,ti,peff_i,zeta_i,Gmem,e_i);
        palisade.add_organelles(plasma_membrane)
        
        #Chloroplast envelop        

        [fi,ti,peff_i,zeta_i,Genv,e_i]=self.inputs.get_chlenvelop()
        chl_envelop = Cell_component(Names.CHLENVELOPE,fi,ti,peff_i,zeta_i,Genv,e_i);
        palisade.add_organelles(chl_envelop)
            
        #Cytosol

        [fi,ti_pal,peff_i,zeta_i,e_i] = self.inputs.get_palisade_cytosol()
        cytosol = Cell_component(Names.CYTOSOL,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        palisade.add_organelles(cytosol)
        
        # Stroma
        [fi,ti_pal,peff_i,zeta_i,e_i]=self.inputs.get_palisade_stroma()
        stroma = Cell_component(Names.STROMA,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        palisade.add_organelles(stroma)
        return palisade
        
    
    def add_spongy_organells(self,spongy):
        physical_constants = Physical_constants()
        [T,P,R,H,DCO2_water] = physical_constants.get_phyical_constants()
        
        """ Add cell components of palisade cells
            plasmalemma and chl. env are already added to palisade
        """
        # cell wall
        [fi,ti_pal,peff_i,zeta_i,e_i] = self.inputs.get_spongy_wall()
        cell_wall = Cell_component(Names.CELLWALL,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        spongy.add_organelles(cell_wall)
                   
        #Cytosol

        [fi,ti_pal,peff_i,zeta_i,e_i] = self.inputs.get_spongy_cytosol()
        cytosol = Cell_component(Names.CYTOSOL,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        spongy.add_organelles(cytosol)
        
        # Stroma
        [fi,ti_pal,peff_i,zeta_i,e_i]=self.inputs.get_spongy_stroma()
        stroma = Cell_component(Names.STROMA,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        spongy.add_organelles(stroma)
        return spongy
    
    def make_cells(self,replicate):
        physical_constants = Physical_constants()
        [T,P,R,H,DCO2_water] = physical_constants.get_phyical_constants()
        
        assumed_values = self.inputs.assumed_values
        [fi,peff_i,zeta_i,e_i]=assumed_values.get_cell_wall()
        
        df = self.inputs.get_leaf_anatomy(replicate)
        pali = df[df['Tissue']=='Palisade']
        spon = df[df['Tissue']=='Spongy']

        Lm_L_tissue_pali = pali[['Lm L']].values[0]
        Lm_L_tissue_spongy = spon[['Lm L']].values[0]
        Lc_Lm_tissue_pali = pali[['Lc Lm']].values[0]
        Lc_Lm_tissue_spongy = spon[['Lc Lm']].values[0]

        gamma_tissue_pali = df[['gamma_tissue_pali']].values[0]
        gamma_tissue_spongy = df[['gamma_tissue_spongy']].values[0]
        f_pal = df[['f_pal']].values[0]

        """ Calculate sub-resistances """

        palisade = Cell(Names.PALISADE,Lm_L_tissue_pali,Lc_Lm_tissue_pali,gamma_tissue_pali,f_pal)

        spongy = Cell(Names.SPONGY,Lm_L_tissue_spongy,Lc_Lm_tissue_spongy,gamma_tissue_spongy,f_pal)


        [fi,ti_pal,peff_i,zeta_i,e_i] = self.inputs.get_palisade_wall(replicate)
        cell_wall = Cell_component(Names.CELLWALL,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        palisade.add_organelles(cell_wall)
        

        [fi,ti_spo,peff_i,zeta_i,e_i] = self.inputs.get_spongy_wall(replicate)
        cell_wall = Cell_component(Names.CELLWALL,fi,ti_spo,peff_i,zeta_i,DCO2_water,e_i);
        spongy.add_organelles(cell_wall)

        # plasma membrane
        [fi,ti,peff_i,zeta_i,Gmem,e_i]=self.inputs.get_plasmallema()
        plasma_membrane = Cell_component(Names.PLASMALLEMA,fi,ti,peff_i,zeta_i,Gmem,e_i);
        palisade.add_organelles(plasma_membrane)
        spongy.add_organelles(plasma_membrane)
    
        #Chloroplast envelop        

        [fi,ti,peff_i,zeta_i,Genv,e_i]=self.inputs.get_chlenvelop()
        chl_envelop = Cell_component(Names.CHLENVELOPE,fi,ti,peff_i,zeta_i,Genv,e_i);
        palisade.add_organelles(chl_envelop)
        spongy.add_organelles(chl_envelop)
            
        #Cytosol

        [fi,ti_pal,peff_i,zeta_i,e_i] = self.inputs.get_palisade_cytosol(replicate)
        cytosol = Cell_component(Names.CYTOSOL,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        palisade.add_organelles(cytosol)

        [fi,ti_spo,peff_i,zeta_i,e_i] = self.inputs.get_spongy_cytosol(replicate)
        cytosol = Cell_component(Names.CYTOSOL,fi,ti_spo,peff_i,zeta_i,DCO2_water,e_i);
        spongy.add_organelles(cytosol)
            
        # Stroma
        [fi,ti_pal,peff_i,zeta_i,e_i]=self.inputs.get_palisade_stroma(replicate)
        stroma = Cell_component(Names.STROMA,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i);
        palisade.add_organelles(stroma)
            
        [fi,ti_spo,peff_i,zeta_i,e_i]=self.inputs.get_spongy_stroma(replicate)
        stroma = Cell_component(Names.STROMA,fi,ti_spo,peff_i,zeta_i,DCO2_water,e_i);
        spongy.add_organelles(stroma)
        
        return [palisade,spongy]
    
    
    def update_cells(self):
        [palisade_old,spongy_old] = self.make_cells()
        if self.cell_component.get_name()==Names.PALISADE:
            palisade = self.cell_component
            palisade = self.add_palisade_organells(palisade)
            spongy=spongy_old                    
        elif self.cell_component.get_name()==Names.SPONGY:
            spongy = self.cell_component
            spongy = self.add_spongy_organells(spongy)
            palisade=palisade_old            
        # cell wall
        elif self.cell_component.get_name() in [Names.CELLWALL,Names.PLASMALLEMA,
                                                  Names.CHLENVELOPE,Names.CYTOSOL,Names.STROMA]:
            if self.types==Names.PALISADE:
                palisade_old.remove_organelles(self.cell_component)
                palisade=palisade_old
                palisade.add_organelles(self.cell_component)
                spongy=spongy_old
            else:
                spongy_old.remove_organelles(self.cell_component)
                spongy=spongy_old
                spongy.add_organelles(self.cell_component)
                palisade=palisade_old
        return [palisade,spongy]