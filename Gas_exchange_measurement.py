"""
Created on Wed Jan 20 18:28:34 2021

@author: Moges
Measurement data of Ci, Iinc
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

PATH = (r'\\WURNET.NL\Homes\retta001\My Documents\Project\2022\1-D Model\CO2_Resistances\Resistance_to_CO2_diffusion_C3_leaf\\')


class Gas_exchange_measurement:
    data = pd.read_excel ('Gas_Exchange_data.xlsx') 
    FORMAT = ['Replicate','Species','Treatment','Measurement_type','Oxygen_level',\
              'Net_CO2_assimilation_rate', 'Intercellular_CO2_concentration', \
                  'PhiPS2','Irradiance','Stomatal_conductance_for_CO2','CO2R']
    df_selected = data[FORMAT]
    A_CI = df_selected.query('Measurement_type=="A_CI_curve"')
    A_I = df_selected.query('Measurement_type=="A_I_curve"')

    def __init__(self,O2,species,treatment):
            self.O2=O2
            self.species=species        
            self.treatment=treatment
            
     
    def set_O2(self,O2):
        self.O2=O2


    def set_species(self,species):
        self.species=species
    
    
    def set_treatment(self,treatment):
        self.treatment=treatment        
         
        
    def get_O2(self):
        return self.O2
      
        
    def get_A_Ci(self):
        return self.A_CI
    
    
    def get_A_I(self):
        return self.A_I


    def get_species(self):
        return self.species
    
    
    def get_treatment(self):
        return self.treatment

    def make_avareage_data(self):
        columns = ['Species','Treatment','Response','Ci','A','Iinc','PhiPS2','Std.dev A','gs','Std.dev gs','Std.dev PhiPS2']
        species = ['BNigra','HIncana']
        treatments = ['HL','LL']
        ave_gas_Exchange_data = pd.DataFrame([])
        df = pd.DataFrame([],columns=columns )
        O = 0.21
        
        # Make average data of ACI and AI
        for plant in species:
                for treatment in treatments:
                    gas_exch_measurement = Gas_exchange_measurement(O,plant,treatment)
                    df_ave = gas_exch_measurement.average_A_CI()
                    
                    I_ave_ci = df_ave['Irradiance'].values
                    Ci_ave_ci = df_ave['Intercellular_CO2_concentration'].values
                    A_ave_ci = df_ave['Net_CO2_assimilation_rate'].values
                    gs_ave_ci = df_ave['Stomatal_conductance_for_CO2'].values
                    phiPS2 = df_ave['PhiPS2'].values
                    A_std = df_ave['Photo_err'].values
                    gs_std = df_ave['gs_err'].values
                    phiPS2_std = df_ave['PhiPS2_err'].values
                    
                    df = pd.DataFrame([],columns=columns )            
                    df['Ci']=Ci_ave_ci; df['A']=A_ave_ci; df['Iinc']=I_ave_ci; 
                    df['Std.dev A']=A_std; df['gs']=gs_ave_ci;df['PhiPS2']=phiPS2;df['Std.dev PhiPS2']=phiPS2_std;
                    df['Std.dev gs']=gs_std;df['Species']=plant; df['Treatment']=treatment; df['Response']='ACI'; 
                    ave_gas_Exchange_data=ave_gas_Exchange_data.append(df)
                    
                    df_ave  = gas_exch_measurement.average_A_I() 
                    I_ave_i = df_ave['Irradiance'].values
                    Ci_ave_i = df_ave['Intercellular_CO2_concentration'].values
                    A_ave_i = df_ave['Net_CO2_assimilation_rate'].values
                    gs_ave_i = df_ave['Stomatal_conductance_for_CO2'].values
                    phiPS2 = df_ave['PhiPS2'].values
                    A_std = df_ave['Photo_err'].values
                    gs_std = df_ave['gs_err'].values
                    phiPS2_std = df_ave['PhiPS2_err'].values
                    
                    
                    df = pd.DataFrame([],columns=columns )
                    df['Ci']=Ci_ave_i; df['A']=A_ave_i; df['Iinc']=I_ave_i; 
                    df['Std.dev A']=A_std; df['gs']=gs_ave_i;df['PhiPS2']=phiPS2;df['Std.dev PhiPS2']=phiPS2_std;
                    df['Std.dev gs']=gs_std;df['Species']=plant; df['Treatment']=treatment; df['Response']='AI'; 
                    ave_gas_Exchange_data=ave_gas_Exchange_data.append(df)
        return ave_gas_Exchange_data
                
    ##ave_gas_Exchange_data.to_excel(PATH + 'Ave_Gas_Exchange_data_corr.xlsx', index = False)


    def get_average_values(self,curve):
        if curve == 'ACI':
            df_ave = self.average_A_CI()    
        else :
            df_ave = self.average_A_I()    
        
        return df_ave


    def plot_A_CI(self):
        A_CI_d = self.A_CI[self.A_CI['Oxygen_level']==self.get_O2()]
        A_CI_d = A_CI_d[A_CI_d['Species']==self.get_species()]
        A_CI_d = A_CI_d[A_CI_d['Treatment']==self.get_treatment()]
        replicates = A_CI_d['Replicate'].values
        replicates=np.unique(replicates)
        for replicate in replicates:
            A_CI_r= A_CI_d[A_CI_d['Replicate']==replicate]
            Ci = A_CI_r['Intercellular_CO2_concentration'].values
            A = A_CI_r['Net_CO2_assimilation_rate'].values
            plt.plot(Ci,A,'o',markersize=8)
            plt.xlabel('Intercellular CO$_2$ (µmol mol$^{-1}$)',fontsize=24)
            plt.ylabel('Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
        plt.show()


    def get_AI_data(self):
        A_I_d  = self.A_I[self.A_I['Oxygen_level']==self.get_O2()]
        A_I_d =  A_I_d[A_I_d['Species']==self.get_species()]
        A_I_d = A_I_d[A_I_d['Treatment']==self.get_treatment()]
        return A_I_d


    def get_ACI_data(self):
        A_CI = self.get_A_Ci()
        A_CI_d  = A_CI[A_CI['Oxygen_level']==self.get_O2()]
        A_CI_d =  A_CI_d[A_CI_d['Species']==self.get_species()]
        A_CI_d = A_CI_d[A_CI_d['Treatment']==self.get_treatment()]
        return A_CI_d

            
    def plot_A_I(self):
        A_I_d  = self.A_I[self.A_I['Oxygen_level']==self.get_O2()]
        A_I_d =  A_I_d[A_I_d['Species']==self.get_species()]
        A_I_d = A_I_d[A_I_d['Treatment']==self.get_treatment()]
        replicates = A_I_d['Replicate'].values
        replicates=np.unique(replicates)
        for replicate in replicates:
            A_I_r= A_I_d[A_I_d['Replicate']==replicate]
            I = A_I_r['Irradiance'].values
            A = A_I_r['Net_CO2_assimilation_rate'].values
            plt.plot(I,A,'o',markersize=8)
            plt.xlabel('Irradiance (µmol m$^{-2}$ s$^{-1}$)',fontsize=24)
            plt.ylabel('Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)',fontsize=24)
            plt.title("A-I")
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
        plt.show()


    def average_A_CI(self):
        A_CI_d = self.A_CI[self.A_CI['Oxygen_level']==self.get_O2()]
        A_CI_d = A_CI_d[A_CI_d['Species']==self.get_species()]
        A_CI_d = A_CI_d[A_CI_d['Treatment']==self.get_treatment()]
        replicates = A_CI_d['Replicate'].unique()

        cols = ['Irradiance','Intercellular_CO2_concentration','Net_CO2_assimilation_rate',\
                'PhiPS2','Stomatal_conductance_for_CO2','Photo_err','gs_err','PhiPS2_err']
        df_ave = pd.DataFrame([],columns = cols)
        df_ci = pd.DataFrame([])
        df_A = pd.DataFrame([])
        df_I = pd.DataFrame([])
        df_gs = pd.DataFrame([])
        df_phi = pd.DataFrame([])
        count = 0

        for replicate in replicates:
            A_CI_r= A_CI_d[A_CI_d['Replicate']==replicate]
            Ci = A_CI_r['Intercellular_CO2_concentration'].values
            A = A_CI_r['Net_CO2_assimilation_rate'].values
            I = A_CI_r['Irradiance'].values
            gs = A_CI_r['Stomatal_conductance_for_CO2'].values
            PhiPS2 = A_CI_r['PhiPS2'].values
            df_ci.loc[:,count] = Ci
            df_A.loc[:,count] = A
            df_I.loc[:,count] = I
            df_gs.loc[:,count] = gs
            df_phi.loc[:,count] = PhiPS2

            count+=1
        
        df_ave.loc[:,'Irradiance'] = np.nanmean(df_I,axis=1)
        df_ave.loc[:,'Intercellular_CO2_concentration'] = np.nanmean(df_ci,axis=1)
        df_ave.loc[:,'Net_CO2_assimilation_rate'] = np.nanmean(df_A,axis=1)
        df_ave.loc[:,'Stomatal_conductance_for_CO2'] = np.nanmean(df_gs,axis=1)
        df_ave.loc[:,'PhiPS2'] = np.nanmean(df_phi,axis=1)
        df_ave.loc[:,'Photo_err'] = np.nanstd(df_A,axis=1)
        df_ave.loc[:,'gs_err'] = np.nanstd(df_gs,axis=1)
        df_ave.loc[:,'PhiPS2_err'] = np.nanstd(df_phi,axis=1)
        df_ave = df_ave.sort_values(by=['Intercellular_CO2_concentration'])
        return df_ave
    
    def average_A_I(self):
        A_I_d  = self.A_I[self.A_I['Oxygen_level']==self.get_O2()]
        A_I_d =  A_I_d[A_I_d['Species']==self.get_species()]
        A_I_d = A_I_d[A_I_d['Treatment']==self.get_treatment()]
        replicates = A_I_d['Replicate'].unique()

        df_ci = pd.DataFrame([])
        df_A = pd.DataFrame([])
        df_I = pd.DataFrame([])
        df_gs = pd.DataFrame([])
        df_phi = pd.DataFrame([])
        count = 0
        cols = ['Irradiance','Intercellular_CO2_concentration','Net_CO2_assimilation_rate',\
                'PhiPS2','Stomatal_conductance_for_CO2','Photo_err','gs_err','PhiPS2_err']
        df_ave = pd.DataFrame([],columns = cols)
        
        for replicate in replicates:
            A_I_r= A_I_d[A_I_d['Replicate']==replicate]
            I = A_I_r['Irradiance'].values
            A = A_I_r['Net_CO2_assimilation_rate'].values
            Ci = A_I_r['Intercellular_CO2_concentration'].values
            gs = A_I_r['Stomatal_conductance_for_CO2'].values
            PhiPS2 = A_I_r['PhiPS2'].values
            df_ci.loc[:,count] = Ci
            df_A.loc[:,count] = A
            df_I.loc[:,count] = I
            df_gs.loc[:,count] = gs
            df_phi.loc[:,count] = PhiPS2
            count+=1
            
        df_ave.loc[:,'Irradiance'] = np.nanmean(df_I,axis=1)
        df_ave.loc[:,'Intercellular_CO2_concentration'] = np.nanmean(df_ci,axis=1)
        df_ave.loc[:,'Net_CO2_assimilation_rate'] = np.nanmean(df_A,axis=1)
        df_ave.loc[:,'Stomatal_conductance_for_CO2'] = np.nanmean(df_gs,axis=1)
        df_ave.loc[:,'PhiPS2'] = np.nanmean(df_phi,axis=1)
        df_ave.loc[:,'Photo_err'] = np.nanstd(df_A,axis=1)
        df_ave.loc[:,'gs_err'] = np.nanstd(df_gs,axis=1)
        df_ave.loc[:,'PhiPS2_err'] = np.nanstd(df_phi,axis=1)        
        df_ave = df_ave.sort_values(by=['Irradiance'])            
        return df_ave
        
    def plot_ave_A_CI(self,df_ave):
        
        Ci_ave_ci = df_ave['Intercellular_CO2_concentration'].values
        A_ave_ci = df_ave['Net_CO2_assimilation_rate'].values
        A_std = df_ave['Photo_err'].values
        plt.errorbar(Ci_ave_ci,A_ave_ci,A_std,fmt='o',markersize=8)
        plt.xlabel('Intercellular CO$_2$ (µmol mol$^{-1}$)',fontsize=24)
        plt.ylabel('Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)',fontsize=24)
        plt.xticks(fontsize=24)
        plt.yticks(fontsize=24)        
        plt.show()



    def plot_ave_A_I(self,df_ave):
        I_ave = df_ave['Irradiance'].values
        A_ave = df_ave['Net_CO2_assimilation_rate'].values
        A_std = df_ave['Photo_err'].values

        plt.errorbar(I_ave,A_ave,A_std,fmt='o')
        plt.xlabel('Irradiance (µmol $m^{-2}$ $s^{-1}$)')
        plt.ylabel('Net photosynthesis (µmol $m^{-2}$ $s^{-1}$)')
        plt.xticks(fontsize=24)
        plt.yticks(fontsize=24)        
        plt.show()


    def plot_ave_gs_I(self,df_ave):
        gs_ave = df_ave['Stomatal_conductance_for_CO2'].values
        gs_std = df_ave['gs_err'].values
        I_ave = df_ave['Irradiance'].values

        plt.errorbar(I_ave,gs_ave,gs_std,fmt='o')
        plt.xlabel('Irradiance (µmol $m^{-2}$ $s^{-1}$)')
        plt.ylabel('Stomatal conductance (mol $m^{-2}$ $s^{-1}$)')
        plt.xticks(fontsize=24)
        plt.yticks(fontsize=24)        
        plt.show()        


    def plot_ave_gs_CI(self,df_ave):
        gs_ave = df_ave['Stomatal_conductance_for_CO2'].values
        gs_std = df_ave['gs_err'].values
        Ci_ave_ci = df_ave['Intercellular_CO2_concentration'].values

        plt.errorbar(Ci_ave_ci,gs_ave,gs_std,fmt='o')
        plt.xlabel('Intercellular $CO_2$ (µmol $mol^{-1}$)')
        plt.ylabel('Stomatal conductance (mol $m^{-2}$ $s^{-1}$)')
        plt.xticks(fontsize=24)
        plt.yticks(fontsize=24)        
        plt.show()
    
 