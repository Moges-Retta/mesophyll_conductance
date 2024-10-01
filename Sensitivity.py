# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 23:20:19 2021

@author: Moges
Sensitivity analysis
"""
from Plant import Plant
from Names import Names
from Data import Data
from Photosynthesis import Photosynthesis
import matplotlib.pyplot as plt
import numpy as np
from Make_mesophyll import Make_mesophyll
from Cell_component import Cell_component
from Cell import Cell
from matplotlib.ticker import MultipleLocator
from Physical_constants import Physical_constants
import pandas as pd
from tabulate import tabulate
# inputs = Data()
physical_constants = Physical_constants()
[T,P,R,H,DCO2_water] = physical_constants.get_phyical_constants()

PATH = (r'\\WURNET.NL\Homes\retta001\My Documents\Project\2021\CO2_Resistances\H.incanaB.nigra\Resistance_to_CO2_diffusion_C3_leaf\\')

class Sensitivity:
    def __init__(self,plant,measurement):
        self.plant=plant
        self.measurement=measurement    
        self.inputs = Data(plant)
        
        
    def calculate_response(self,plant,curve):
        """ Calculate response curve.
        
        Parameters:
        -----------------
        
        arg_1 : Object
        
        Plant object
        
        arg_2 : str
        A-CI or A-Iinc curve
        
        Returns:
        ----------------
        Array of floats
                
        AN_CI_mod : A-CI response
        
        F_CI : Response of photorespiration
        
        f_refix : fraction of total refixation, cell = cell, 
                ias, intercellular air space   
        """
        
        f_refix_ci =[]; AN_CI_mod=[];F_CI=[];f_refix_cell=[];f_refix_ias=[];
    
        if curve == Names.LIGHT:
            df_ave =  self.measurement.average_A_I()
            I_ave_i = df_ave['Irradiance'].values
            Ci_ave_i = df_ave['Intercellular_CO2_concentration'].values
            A_ave_i = df_ave['Net_CO2_assimilation_rate'].values
            gs_ave_i = df_ave['Stomatal_conductance_for_CO2'].values

            i=0
            for Iinc in I_ave_i:
                photosynthesis = Photosynthesis(plant,Iinc,Ci_ave_i[i],self.measurement.get_O2(),gs_ave_i[i])
                [A,F,Cc,rcx] = photosynthesis.calculate_A()
                [f,f_cell,f_ias] = photosynthesis.calculate_f_refix() 
                AN_CI_mod.append(A)
                F_CI.append(F)
                f_refix_ci.append(f)
                f_refix_cell.append(f_refix_cell)
                f_refix_ias.append(f_refix_ias)
                i+=1
        else:
            df_ave =  self.measurement.average_A_CI()

            I_ave_ci = df_ave['Irradiance'].values
            Ci_ave_ci = df_ave['Intercellular_CO2_concentration'].values
            A_ave_ci = df_ave['Net_CO2_assimilation_rate'].values
            gs_ave_ci = df_ave['Stomatal_conductance_for_CO2'].values
            i=0
            for Ci in Ci_ave_ci:
                photosynthesis = Photosynthesis(plant,I_ave_ci[i],Ci,self.measurement.get_O2(),gs_ave_ci[i])
                [A,F,Cc,rcx] = photosynthesis.calculate_A()
                [f,f_cell,f_ias] = photosynthesis.calculate_f_refix() 
                AN_CI_mod.append(A)
                F_CI.append(F)
                f_refix_ci.append(f)
                f_refix_cell.append(f_refix_cell)
                f_refix_ias.append(f_refix_ias)
                i+=1
        return [AN_CI_mod,F_CI,f_refix_ci,f_refix_cell,f_refix_ias]


    def plot_sensitivity(self,AN_CI_mod_sensitivity,AN_I_mod_sensitivity,Ci_ave_ci,I_ave_i,title,A_ave_i,A_ave_ci):

        fig, ax = plt.subplots(1,2,constrained_layout=False)
        ax[0].plot(Ci_ave_ci, AN_CI_mod_sensitivity[:,0],'k-.', label='-25%', 
                linewidth=1.4, markersize=10)
        plt.rcParams["figure.figsize"] = (15,10)
        plt.rcParams.update({'font.size': 16})
        ax[0].plot(Ci_ave_ci, AN_CI_mod_sensitivity[:,1], 'k', label='default.')
        ax[0].plot(Ci_ave_ci, AN_CI_mod_sensitivity[:,2], 'k--.', label='+25%')
        ax[0].plot(Ci_ave_ci, A_ave_ci[:], 'o', label='expt')
        
        ax[0].tick_params(labelsize='medium', width=2)
        ax[0].xaxis.set_minor_locator(MultipleLocator(250))
        ax[0].yaxis.set_minor_locator(MultipleLocator(5))
        ax[0].tick_params(which='major', length=7)
        ax[0].tick_params(which='minor', length=5)
        ax[0].legend(loc='lower right', fontsize='x-large')
        ax[0].set_ylabel("Net photosynthesis (µmol $m^{-2}$ $s^{-1}$)",fontsize=24)
        ax[0].set_xlabel("Intercellular $CO_2$ (µmol $mol^{-1}$)",fontsize=24)
        #plt.savefig("A_CI.tif", dpi=300)
        #plot A-I:

        plt.rcParams["figure.figsize"] = (15,10)
        plt.rcParams.update({'font.size': 16})
        ax[1].plot(I_ave_i, AN_I_mod_sensitivity[:,0],'k-.', label='-25%')
        ax[1].plot(I_ave_i, AN_I_mod_sensitivity[:,1], 'k', label='default')
        ax[1].plot(I_ave_i, AN_I_mod_sensitivity[:,2], 'k--.', label='+25%')
        ax[1].plot(I_ave_i, A_ave_i[:], 'o', label='expt')
        
        ax[1].tick_params(labelsize='medium', width=2)
        # ax[1].legend(loc='lower right',fontsize='x-large')
        ax[1].tick_params(which='major', length=7)
        ax[1].tick_params(which='minor', length=5)
        ax[1].xaxis.set_minor_locator(MultipleLocator(125))
        ax[1].yaxis.set_minor_locator(MultipleLocator(2.5))
        ax[1].set_ylabel("Net photosynthesis (µmol $m^{-2}$ $s^{-1}$)",fontsize=24)
        ax[1].set_xlabel("Irradiance (µmol $m^{-2}$ $s^{-1}$)",fontsize=24)  
        fig.suptitle(title,ha = 'left', va='bottom')
        plt.tight_layout()
        plt.show()


# Calculate the response of AN-CI or AN-I to changes in leaf anatomy
    def calculate_sensitivity(self,item,params,whichParm,Ci_ave_ci,I_ave_ci,curve,name_component,cell_type_code):
        inputs = Data(self.plant)


        if curve ==Names.CO2:
            AN_CI_mod_sensitivity=np.zeros((len(Ci_ave_ci),3))
        else:
            AN_CI_mod_sensitivity=np.zeros((len(I_ave_ci),3))
        j=0
        factors=[0.25,1,1.25]
        old_param=params[whichParm]
        df = pd.DataFrame([],columns=('Tissue','Cell component','parameter','factor','rdiff'))
        for thisParam in [i*params[whichParm] for i in factors]:
            params[whichParm] = thisParam
            if item.get_name()==Names.PALISADE:
                palisade = Cell(item.get_name(),params[0],params[1],params[2],params[3])
                myMesophyll = Make_mesophyll(inputs,palisade,Names.PALISADE)
            elif item.get_name()==Names.SPONGY:
                spongy = Cell(item.get_name(),params[0],params[1],params[2],params[3])
                myMesophyll = Make_mesophyll(inputs,spongy,Names.SPONGY)
            else:
                cell = Cell_component(item.get_name(),params[0],params[1],params[2],params[3],params[4]);
                if cell_type_code==Names(Names.PALISADE).value:
                    myMesophyll = Make_mesophyll(inputs,cell,Names.PALISADE)
                else:
                    myMesophyll = Make_mesophyll(inputs,cell,Names.SPONGY)
            params[whichParm] = old_param   
            [palisade,spongy] = myMesophyll.update_cells()
            admiro = Plant(self.plant.get_name(),self.plant.get_treatment(),[self.plant.get_lambda_mitochondria(),self.plant.get_k()])
            admiro.add_cells(palisade)
            admiro.add_cells(spongy)
            rdiff= admiro.calculate_rdiff()
            rdiff= rdiff*H/10**5
            [AN_CI_mod,F_CI,f_refix_ci,f_refix_cell,f_refix_ias] = self.calculate_response(admiro,curve)
            if curve != Names.CO2:
                 df.loc[j,'factor']=factors[j]
                 df.loc[j,'rdiff']=rdiff[0]
                 df.loc[j,'Cell component']=Names(item.get_name()).value;
                 df.loc[j,'parameter']=name_component;
                 df.loc[j,'Tissue']=cell_type_code
                 # df=df.append([tam])
            AN_CI_mod_sensitivity[:,j]= AN_CI_mod
            j+=1
        if j==3:
            if not df.empty:
                print(tabulate(df, headers='keys', tablefmt='psql',showindex=False))
        return [AN_CI_mod_sensitivity,df]


    def sensitivity_in_response(self,cell_component,params,titles,cell_type_code):
        plant =self.plant.get_name();
        treatment =self.plant.get_treatment();
        
        df_ave= self.measurement.average_A_I()
        I_ave_i = df_ave['Irradiance'].values
        Ci_ave_i = df_ave['Intercellular_CO2_concentration'].values
        A_ave_i = df_ave['Net_CO2_assimilation_rate'].values
        # gs_ave_i = df_ave['Stomatal_conductance_for_CO2'].values
        
        df_ave = self.measurement.average_A_CI()

        I_ave_ci = df_ave['Irradiance'].values
        Ci_ave_ci = df_ave['Intercellular_CO2_concentration'].values
        A_ave_ci = df_ave['Net_CO2_assimilation_rate'].values
        # gs_ave_ci = df_ave['Stomatal_conductance_for_CO2'].values
        i=0
        df_tot=pd.DataFrame();
        for param in params:
            [AN_CI_mod_sensitivity,df] = self.calculate_sensitivity(cell_component,params,i,Ci_ave_ci,I_ave_ci,Names.CO2,titles[i],cell_type_code)
            [AN_I_mod_sensitivity,df] = self.calculate_sensitivity(cell_component,params,i,Ci_ave_i,I_ave_i,Names.LIGHT,titles[i],cell_type_code)
            self.plot_sensitivity(AN_CI_mod_sensitivity,AN_I_mod_sensitivity,Ci_ave_ci,I_ave_i,titles[i],A_ave_i,A_ave_ci)
            i+=1
            df_tot=df_tot.append(df)
        df_tot.to_excel(PATH + 'Sensitivity_'+plant+'_'+treatment+ '_'+
                cell_type_code+' '+Names(cell_component.get_name()).value+'.xlsx', index = False)
# Cell wall
    def of_cell_wall(self,cell_type_code):
        inputs = Data(self.plant)
        
        if cell_type_code == Names.SPONGY.value:
            [fi,ti_pal,peff_i,zeta_i] = inputs.get_spongy_wall() 
        else:
            [fi,ti_pal,peff_i,zeta_i] = inputs.get_palisade_wall()

        cell_wall = Cell_component(Names.CELLWALL,fi,ti_pal,peff_i,zeta_i,DCO2_water);
        params = [fi,ti_pal,peff_i,zeta_i,DCO2_water]
        titles = ("fi","ti","peff","zeta","DCO2_water")    
        self.sensitivity_in_response(cell_wall,params,titles,cell_type_code)


# stroma
    def of_stroma(self,cell_type_code): 
        inputs = Data(self.plant)
        
        if cell_type_code == Names.SPONGY.value:
            [fi,ti_pal,peff_i,zeta_i] = inputs.get_spongy_stroma()  
        else:
            [fi,ti_pal,peff_i,zeta_i] = inputs.get_palisade_stroma()
            
        stroma = Cell_component(Names.STROMA,fi,ti_pal,peff_i,zeta_i,DCO2_water);
        params = [fi,ti_pal,peff_i,zeta_i,DCO2_water]
        titles = ("fi","ti","peff","zeta","DCO2_water")
        self.sensitivity_in_response(stroma,params,titles,cell_type_code)


    # Chloroplast envelop
    def of_chloroplast_envelop(self,cell_type_code): 
        
        inputs = Data(self.plant)
        if cell_type_code == Names.SPONGY.value:
            [fi,ti_pal,peff_i,zeta_i,Genv] = inputs.get_chlenvelop()  
        else:
            [fi,ti_pal,peff_i,zeta_i,Genv] = inputs.get_chlenvelop()
            
        chloroplast = Cell_component(Names.CHLENVELOPE,fi,ti_pal,peff_i,zeta_i,Genv);
        params = [fi,ti_pal,peff_i,zeta_i,Genv]
        titles = ("fi","ti","peff","zeta","Genv")
        self.sensitivity_in_response(chloroplast,params,titles,cell_type_code)
    
    
    # Cytosol
    def of_cytosol(self,cell_type_code):
        inputs = Data(self.plant)
        if cell_type_code == Names.SPONGY.value:
            [fi,ti_pal,peff_i,zeta_i] = inputs.get_spongy_cytosol()  
        else:
            [fi,ti_pal,peff_i,zeta_i] = inputs.get_palisade_cytosol()  
            
        cytosol = Cell_component(Names.CYTOSOL,fi,ti_pal,peff_i,zeta_i,DCO2_water);
        params = [fi,ti_pal,peff_i,zeta_i,DCO2_water]
        titles = ("fi","ti","peff","zeta","DCO2_water")
        self.sensitivity_in_response(cytosol,params,titles,cell_type_code)
        
        
    ## palisade
    def of_palisade(self,cell_type_code): 
        inputs = Data(self.plant)
        [Lm_L_tissue_pali,Lm_L_tissue_spongy,Lc_Lm_tissue_pali,Lc_Lm_tissue_spongy,
         gamma_tissue_pali,gamma_tissue_spongy,f_pal] = inputs.get_leaf_anatomy()
        palisade = Cell(Names.PALISADE,Lm_L_tissue_pali,Lc_Lm_tissue_pali,gamma_tissue_pali,f_pal)
        params = [Lm_L_tissue_pali,Lc_Lm_tissue_pali,gamma_tissue_pali,f_pal]
        titles = ("Lm_L_tissue_pali","Lc_Lm_tissue_pali","gamma_tissue_pali","f_pal")
        self.sensitivity_in_response(palisade,params,titles,cell_type_code)
        
        
    ## spongy
    def of_spongy(self,cell_type_code):  
        inputs = Data(self.plant)
        [Lm_L_tissue_pali,Lm_L_tissue_spongy,Lc_Lm_tissue_pali,Lc_Lm_tissue_spongy,gamma_tissue_pali,gamma_tissue_spongy,f_pal] = inputs.get_leaf_anatomy()
        spongy = Cell(Names.SPONGY,Lm_L_tissue_spongy,Lc_Lm_tissue_spongy,gamma_tissue_spongy,f_pal)
        params = [Lm_L_tissue_spongy,Lc_Lm_tissue_spongy,gamma_tissue_spongy,f_pal]
        titles = ("Lm_L_tissue_spongy","Lc_Lm_tissue_spongy","gamma_tissue_spongy","f_pal")
        self.sensitivity_in_response(spongy,params,titles,cell_type_code)