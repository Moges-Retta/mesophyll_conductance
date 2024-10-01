# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:54:53 2021
@author: Moges Retta

""" 

from Names import Names
from Data import Data
import pandas as pd
from Physical_constants import Physical_constants
import numpy as np
from tabulate import tabulate
from Assumed_values import Assumed_values

class Plant:
    """
    Model cultivars as a plant with characteristics:
        An : Photosynthesis rate : µmol/m2/s
        rdiff : Resistance to diffusion 
        w : ratio of rchl to  rdiff
    
    Calculate resistance rdiff
         based on Berghuijs et al. 2015, PS,10.1016/j.plantsci.2015.06.022              
    """     

    def __init__(self, name,treatment,mitochondria_location):
        self.name=name
        self.treatment=treatment
        self.Cells = []
        self.lambda_mitochondria = mitochondria_location[0]
        self.k = mitochondria_location[1]
        
        
    def add_cells(self,cell):
            self.Cells.append(cell)
            
    def remove_cells(self,cell):
            self.Cells.remove(cell)
            
    def set_name(self,name):
        if isinstance(name, str):
            self.name = name
   
         
    def set_lambda_mitochondria(self,lambda_mitochondria):
            self.lambda_mitochondria = lambda_mitochondria


    def set_k(self,k):
            self.k = k


    def set_treatment(self,treatment):
            self.treatment = treatment

            
    def get_treatment(self):
        return self.treatment
    
    
    def get_name(self):
        return self.name
    
    
    def get_eTransportConstants(self):
        inputs = Data(self)             
        return  inputs.get_electron()
    
    
    def get_Jmax(self):
        inputs = Data(self)
        return  inputs.get_electron()[0]
    
    
    def get_k2LL(self):
        inputs = Data(self)
        return  inputs.get_electron()[2]
    
    
    def get_teta(self):
        inputs = Data(self)
        return  inputs.get_electron()[1]
    
        
     # [Vcmax,KmC,KmO,Rd,Sco,Tp]      
    def get_rubiscoConstants(self):
        inputs = Data(self)
        return  inputs.get_rubisco_kinetics()
    
    
    def get_Vcmax(self):
        inputs = Data(self)
        return  inputs.get_rubisco_kinetics()[0]
    
    
    def get_KmC(self):
        inputs = Data(self)
        return  inputs.get_rubisco_kinetics()[1]
    
    
    def get_KmO(self):
        inputs = Data(self)
        return  inputs.get_rubisco_kinetics()[2]
    
    
    def get_Rd(self):
        inputs = Data(self)
        return  inputs.get_rubisco_kinetics()[3] 
    
    
    def get_Sco(self):
        inputs = Data(self)
        return  inputs.get_rubisco_kinetics()[4]
    
    
    def get_Tp(self):
        inputs = Data(self)
        return  inputs.get_rubisco_kinetics()[5]
    
    
    def get_cells(self):
            return  self.Cells   
        
        
    def get_k(self):
            return  self.k   
        
        
    def get_lambda_mitochondria(self):
            return  self.lambda_mitochondria   
        
        
    def get_alpha(self):
            return  self.lambda_mitochondria*self.k  
     
        
    def get_porosity(self,replicate):
        inputs = Data(self)
        return  inputs.get_leaf_porosity(replicate)
    
    
    def get_mesophyll_thickness(self,replicate):
        inputs = Data(self)
        return  inputs.get_mesophyll_thickness(replicate)
    
    def get_sigma(self):
        inputs = Data(self)
        return  inputs.get_sigma()
    
    def calculate_rias(self,replicate):
            """ Calculate intercellular airspace.
            
            Parameters:
            -----------------
            
            arg_1 : Object
            Plant object
            
            Returns:
            ----------------
            Float
            
            rias : 
            Resistance for CO2 transport of the intercellular airspace
            
            m2 leaf s bar CO2 mol−1 CO2
            
        
            """
             # = 4 # amphistomatous leaves has effective length of 1/4th
            dCO2_gas = Physical_constants.DCO2_gas
            turtosity = Assumed_values.turtosity
            porosity=self.get_porosity(replicate)
            l_ias=self.get_mesophyll_thickness(replicate)*Assumed_values.path_lengthening/Assumed_values.factor
            return l_ias*turtosity/(dCO2_gas*porosity)
                     
                       
           
    def calculate_rliq(self):
        """ Calculate rdiff.
        
        Parameters:
        -----------------
        
        arg_1 : Object
        Plant object
        
        Returns:
        ----------------
        Float
        
        rdiff : 
        Total resistance for CO2 transport of the physical barriers 
        in the mesophyll, 
        m2 leaf s bar CO2 mol−1 CO2
        
        w : Ratio of rchl to rdiff
    
        """
        
        rdiff = 0;
        rc=0;
        for cell in self.Cells:
            for component in cell.get_cell_components():
                r = component.calculate_component_resistance()
                if cell.get_name()==Names.PALISADE:    
                   # if component.get_name() in [Names.STROMA,Names.CELLWALL,Names.CYTOSOL]:
                    r *= cell.get_f_pal()
                else:
                   # if component.get_name() in [Names.STROMA,Names.CELLWALL,Names.CYTOSOL]:
                    r *= (1 - cell.get_f_pal())
                if component.get_name() in [Names.STROMA,Names.CHLENVELOPE]:
                    r /= self.calculate_Sc_S()/Assumed_values.connectivity;
                else:
                    r /= self.calculate_Sm_S()/Assumed_values.connectivity;
                    
                if component.get_name() in [Names.STROMA,Names.CHLENVELOPE,Names.CYTOSOL]:
                    if component.get_name()==Names.CYTOSOL:
                        rc+=r/2
                    else:
                        rc+=r                   
                rdiff+=r
        return [rdiff, rc/rdiff]
    

    def calculate_rdiff(self):
                """ Calculate rdiff.
                
                Parameters:
                -----------------
                
                arg_1 : Object
                Plant object
                
                Returns:
                ----------------
                Float
                
                rias : 
                Total resistance for CO2 transport of the physical barriers 
                        in the leaf                
                m2 leaf s bar CO2 mol−1 CO2
                
                """
                [rliq,w] = self.calculate_rliq();
                return rliq*Physical_constants.R*Physical_constants.T/Physical_constants.H + self.calculate_rias(4);  #replicate 4 is average value
            
    def calculate_r_individual(self,replicate):
        """ Calculate resistance of individual cell components.
        
        Parameters:
        -----------------
        
        arg_1 : Object
        Plant object
        
        Returns:
        ----------------
        data frame
        
        r : 
        Resistance for CO2 transport of the physical barriers 
        in the mesophyll, 
        m2 leaf s bar CO2 mol−1 CO2
           
        """
        df = pd.DataFrame([],columns=['Replicate','Cell type','Cell component','Resistance','% tot','w']);
        j=0
        rc = 0
        rdiff = 0;
        
        for cell in self.Cells:
            for component in cell.get_cell_components():
                r = component.calculate_component_resistance()
                
                # if component.get_name() in [Names.STROMA, Names.CELLWALL, Names.CYTOSOL]:
                if cell.get_name() == Names.PALISADE:
                    r *= cell.get_f_pal()
                else:
                    r *= (1 - cell.get_f_pal())
                    
                if component.get_name() in [Names.STROMA, Names.CHLENVELOPE]:
                    r /= self.calculate_Sc_S()/Assumed_values.connectivity;
                else:
                    r /= self.calculate_Sm_S()/Assumed_values.connectivity;
    
                if component.get_name() in [Names.STROMA,Names.CHLENVELOPE,Names.CYTOSOL]:
                    if component.get_name()==Names.CYTOSOL:
                        rc+=r/2
                    else:
                        rc+=r    
                rdiff+=r
        
                df.loc[j,'Cell component']=Names(component.get_name()).name
                df.loc[j,'Resistance']=np.round(r[0]*Physical_constants.H/10**5,3)
                df.loc[j,'Cell type']=Names(cell.get_name()).name
                j+=1
                
        w = rc/rdiff
        df.loc[:,'w'] = w[0]
            
        tot_val=df['Resistance'].values
        tot=tot_val.sum()
        df.loc[:,'% tot']=(tot_val/tot)*100
        df.loc[:,'% tot']=np.round(df.loc[:,'% tot'].values.astype(float))
        df.loc[:,'Replicate']=replicate
 
        # df=df.sort_values(by='% tot', ascending=False)
        
        # print(tabulate(df, headers='keys', tablefmt='psql',showindex=False))
        
        

        return df
    
    def calculate_Sc_S(self):
        """ Calculate Sc/S. based on Eq. 4-5
        
        Parameters:
        -----------------
        
        arg_1 : Object
        Plant object
        
        Returns:
        ----------------
        Float
        
        Fraction of the exposed chloroplast surface area of the 
        palisade parenchyma and the spongy parenchyma relative to 
        leaf surface area at one one side of the leaf, m/m
    
        """
        Sc_S = 0
        for cell in self.Cells:
            Sc_S += cell.calculate_Sc_S();
        return Sc_S


    def calculate_Sm_S(self):
        """ Calculate Sm/S.
        
        Parameters:
        -----------------
        
        arg_1 : Object
        Plant object
        
        Returns:
        ----------------
        Float
        
        Fraction of the exposed mesophyll surface area of the
        palisade parenchyma and the spongy parenchyma relative to 
        leaf surface area at one side of the leaf, m/m
    
        """
        Sm_S = 0
        for cell in self.Cells:
            Sm_S += cell.calculate_Sm_S();
        return Sm_S