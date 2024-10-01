# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 15:48:43 2021

@author: Moges Retta
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:34:12 2021

@author: Moges Retta
Calculate resistances to CO2 transport by leaf microstructure
Analyse sensitivity of rate of photosynthesis to changes in leaf anatomy
based on Berghuijs et al. 2015, PS,10.1016/j.plantsci.2015.06.022 


        Variables
 ----------------------------------------------------------------------
DCO2_i          # Diffusion coefficient of CO2 in component i, 
                  m2/s

DCO2,water      # Diffusion coefficient of CO2 in water, 
                  m2/s

fi              # fraction of the diffusive path length of component and its 
                  thickness, 
                  m/m

fpal            # raction of the exposed mesophyll surface area that belongs 
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
from Plant import Plant
from Names import Names
from Data import Data
from Gas_exchange_measurement import Gas_exchange_measurement
from Photosynthesis import Photosynthesis
import matplotlib.pyplot as plt
import numpy as np
from Make_mesophyll import Make_mesophyll
from Cell_component import Cell_component
from Sensitivity import Sensitivity
from Physical_constants import Physical_constants

def rdiff_value(plant):
    """ Constants"""
    physical_constants = Physical_constants()
    [T,P,R,H,DCO2_water] = physical_constants.get_phyical_constants()
    
    parameters = Data(plant)
    
    [fi,ti_pal,peff_i,zeta_i,e_i] = parameters.get_palisade_wall()
    
    cell_wall = Cell_component(Names.CELLWALL,fi,ti_pal,peff_i,zeta_i,DCO2_water)
    myMesophyll = Make_mesophyll(parameters,cell_wall,"")
    [palisade,spongy] = myMesophyll.make_cells()
    
    """ calculate rdiff, An and w  """
    plant.add_cells(palisade)
    [rdiff_pali,w] = plant.calculate_rdiff()
    rdiff_pali=rdiff_pali*H/10**5
    Sc_S_pali = plant.calculate_Sc_S()
    Sm_S_pali = plant.calculate_Sm_S()
    
    plant.add_cells(spongy)
    [rdiff,w] = plant.calculate_rdiff()
    rdiff= rdiff*H/10**5
    Sc_S_tissue = plant.calculate_Sc_S()
    Sm_S_tissue = plant.calculate_Sm_S()
    
    Sc_S_spongy = Sc_S_tissue - Sc_S_pali
    rdiff_spongy = rdiff-rdiff_pali
    Sm_S_spongy = Sm_S_tissue - Sm_S_pali
    print(plant.calculate_r_individual())

    print("Tis. Sm/s Sc/S")
    print("pal" + " " + str(round(Sm_S_pali[0],2)) + " " + str(round(Sc_S_pali[0],2)))
    print("spo" + " "+ str(round(Sm_S_spongy[0],2))+ " " + str(round(Sc_S_spongy[0],2)))
    print("mes" + " " + str(round(Sm_S_tissue[0],2))+ " "+ str(round(Sc_S_tissue[0],2)))

    print("rdiff:" + " "+ str(round(rdiff[0],3)))
    print("w:" + " "+ str(round(w[0],3)))

""" Photosynthesis, photorespiration and refixation"""

def calculate_response(curve,plant,df_ave):
    I_ave = df_ave['Irradiance'].values
    Ci_ave = df_ave['Intercellular_CO2_concentration'].values
    gs_ave = df_ave['Net_CO2_assimilation_rate'].values
    O = 210

    f_refix_ci =[]; F_CI=[];f_refix_cell=[];f_refix_ias=[];
    AN_CI_mod=[];
    
    I =  [Iinc for Iinc in I_ave] 
    CI =  [Ci for Ci in Ci_ave]   
    if curve == Names.LIGHT:
        photosynthesis_obj =  [Photosynthesis(plant,Iinc,Ci_ave[I.index(Iinc)],O,gs_ave[I.index(Iinc)]) for Iinc in I_ave] 
    else:
        photosynthesis_obj =  [Photosynthesis(plant,I_ave[CI.index(Ci)],Ci,O,gs_ave[CI.index(Ci)]) for Ci in Ci_ave] 

    x = [photosynthesis.calculate_A() for photosynthesis in photosynthesis_obj] 
    f = [photosynthesis.calculate_f_refix() for photosynthesis in photosynthesis_obj]
    x = np.asarray(x)
    AN_CI_mod.append(x[:,0])
    F_CI.append(x[:,1])
    f=np.asarray(f)
    f_refix_ci.append(f[:,0])
    f_refix_cell.append(f[:,1])
    f_refix_ias.append(f[:,2])
    return [AN_CI_mod,F_CI,f_refix_ci,f_refix_cell,f_refix_ias]


# calculate_refix
def calculate_refix(plant,df_ave):
    I_ave = df_ave['Irradiance'].values
    Ci_ave = df_ave['Intercellular_CO2_concentration'].values
    gs_ave = df_ave['Net_CO2_assimilation_rate'].values
    
    I = I_ave[4]
    Ci = Ci_ave[4]
    gs= gs_ave[4]
    O = 210
    photosynthesis = Photosynthesis(plant,I,Ci,O,gs)
    x = photosynthesis.calculate_A()
    f = photosynthesis.calculate_f_refix()
    x = np.asarray(x)
    return [x,f]

def plot_response(curve,df_ave,y_mod):
    y=df_ave['Net_CO2_assimilation_rate'].values
    fig, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = (10,8)
    plt.rcParams.update({'font.size': 16})
    plt.ylabel("Net photosynthesis (µmol $m^{-2}$ $s^{-1}$)",fontsize=24)
    if curve==Names.CO2:
        x=df_ave['Intercellular_CO2_concentration'].values
        ax.plot(x, y_mod,'k', label='Mod.', 
            linewidth=1.4, markersize=10)
        ax.plot(x, y, 'ko', label='Expt.')
        plt.xlabel("Intercellular $CO_2$ (µmol $mol^{-1}$)",fontsize=24)  
    else:
        x = df_ave['Irradiance'].values
        ax.plot(x, y_mod,'k', label='Mod.')
        ax.plot(x, y, 'ko', label='Expt.')
        plt.xlabel("Irradiance (µmol $mol^{-1}$)",fontsize=24)  
    ax.tick_params(labelsize='medium', width=2)
    ax.legend(loc='lower right', fontsize='x-large')

def plot_refixation(df_ave_ci,df_ave_i,f_refix_ci,f_refix_cell_ci,f_refix_ias_ci,f_refix_i,f_refix_cell_i,f_refix_ias_i):
    Ci_ave_ci = df_ave_ci['Intercellular_CO2_concentration'].values
    I_ave_i = df_ave_i['Irradiance'].values
    
    plt.rcParams["figure.figsize"] = (10,8)
    plt.rcParams.update({'font.size': 16})
    fig, ax = plt.subplots(1,2, sharey=True)
    ax[0].plot(Ci_ave_ci.transpose(),np.asarray(f_refix_ci)[0],'ko', label='$f_{refix}$',markersize=12)
    ax[0].plot((Ci_ave_ci.transpose()),np.asarray(f_refix_cell_ci)[0],'k^', label='$f_{refix,cell}$',fillstyle='none',markersize=12)
    ax[0].plot((Ci_ave_ci.transpose()),np.asarray(f_refix_ias_ci)[0],'ks', label='$f_{refix,ias}$',fillstyle='none',markersize=12)
    ax[0].set_ylabel("Fraction of refixation (-)")
    ax[0].set_xlabel("Intercellular $CO_2$ (µmol $mol^{-1}$)")  
    ax[0].set_ylim(top=0.5)

    ax[1].plot( I_ave_i,np.asarray(f_refix_i)[0],'ko', label='$f_{refix}$',markersize=12)
    ax[1].plot( I_ave_i,np.asarray(f_refix_cell_i)[0],'k^', label='$f_{refix,cell}$',fillstyle='none',markersize=12)
    ax[1].plot( I_ave_i,np.asarray(f_refix_ias_i)[0],'ks', label='$f_{refix,ias}$',fillstyle='none',markersize=12)
    ax[1].set_xlabel("Irradiance (µmol $mol^{-1}$)")  
    ax[0].legend(loc='upper right')
    # ax[1].set_ylim(top=0.5)
    plt.tight_layout()

# calculate rdiff
""" Plant and Mitochondria chloroplast arrangement"""
gap_effect = Names.GAP_NOCHANGE
mitochondria_location=Names.INNER_CYTOSOL
chloroplast_mitochondria_arrangement = [mitochondria_location.value,gap_effect.value]
    
treatment = 'HL'
plant = Plant(Names.HIncana.value,treatment,chloroplast_mitochondria_arrangement)
rdiff_value(plant)

    
# Photosynthesis and refixation responses
O = 0.21
# gas_exch_measurement = Gas_exchange_measurement(O)
gas_exch_measurement = Gas_exchange_measurement(O,plant.get_name(),treatment)    
df_ave_ci = gas_exch_measurement.average_A_CI()

[AN_CI_mod,F_CI,f_refix_ci,f_refix_cell_ci,f_refix_ias_ci] = calculate_response(Names.CO2, plant,df_ave_ci)
[A,F]=calculate_refix(plant,df_ave_ci)
plot_response(Names.CO2,df_ave_ci,np.asarray(AN_CI_mod).flatten())

df_ave_i= gas_exch_measurement.average_A_I()
[AN_I_mod,F_I,f_refix_i,f_refix_cell_i,f_refix_ias_i] = calculate_response(Names.LIGHT,plant,df_ave_i)

plot_response(Names.LIGHT,df_ave_i,np.asarray(AN_I_mod).flatten())
plot_refixation(df_ave_ci,df_ave_i,f_refix_ci,f_refix_cell_ci,f_refix_ias_ci,f_refix_i,f_refix_cell_i,f_refix_ias_i) 

np.asarray(f_refix_ci).flatten()[5]

print("fraction of refixed CO2, 400 µmol/mol:" + " "+ str(np.round(np.asarray(f_refix_ci).flatten()[5],3)))
    
#""" Sensitivity """
#sensitivity = Sensitivity(incana,gas_exch_measurement)
#sensitivity.of_cell_wall()
#sensitivity.of_chloroplast_envelop()
#sensitivity.of_cytosol()
#sensitivity.of_stroma()
#sensitivity.of_palisade()
#sensitivity.of_spongy()

# Analyse sensitivity
import pandas as pd
import seaborn as sns

PATH = (r'\\WURNET.NL\Homes\retta001\My Documents\Project\2021\CO2_Resistances\H.incanaB.nigra\Resistance_to_CO2_diffusion_C3_leaf\\')

cell_wall = pd.read_excel ('Sensitivity_H.Incana_HL_Palisade Cell wall.xlsx') 
chl_env = pd.read_excel ('Sensitivity_H.Incana_HL_Palisade Chl. envelope.xlsx') 
cyto = pd.read_excel ('Sensitivity_H.Incana_HL_Palisade Cytosol.xlsx') 
stroma = pd.read_excel ('Sensitivity_H.Incana_HL_Palisade Stroma.xlsx')
palisade = pd.read_excel ('Sensitivity_H.Incana_HL_Palisade Palisade.xlsx')
spongy = pd.read_excel ('Sensitivity_H.Incana_HL_Palisade Spongy.xlsx')

df=pd.DataFrame([])
df = df.append(cell_wall); df = df.append(chl_env);df = df.append(cyto);
df = df.append(stroma)
df.to_excel(PATH+'sensitivity_Hi_Hl_palisade_all.xlsx')

                     
ax1=sns.scatterplot(x="parameter", y="rdiff", hue='factor', size='factor',
                    sizes=(50, 100), legend="full",data=df,palette="deep")
plt.setp(ax1.spines.values(), linewidth=1)

plt.xlabel(" ", fontsize= 24)
plt.ylabel("r$_{diff}$ (bar m${^2}$ s mol$^{-1}$)", fontsize= 24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)


df=pd.DataFrame([])
df = df.append(palisade); df = df.append(spongy);
                     
ax1=sns.scatterplot(x="parameter", y="rdiff", hue='factor', size='factor',
                    sizes=(50, 100), legend="full",data=df,palette="deep")
plt.setp(ax1.spines.values(), linewidth=1)

plt.xlabel(" ", fontsize= 24,rotation='vertical')
plt.ylabel("r$_{diff}$ (bar m${^2}$ s mol$^{-1}$)", fontsize= 24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.setp(ax1.get_xticklabels(), rotation=45)
ax1.set_xticklabels(['Lc_pali','Lm_pali','CCF_pali','f_pal','Lc_spo','Lm_spo','CCF_spo'])



HiHL = pd.read_excel ('Refixation_H.Incana_HL_.xlsx') 
HiLL = pd.read_excel ('Refixation_H.Incana_LL_.xlsx') 
BnHL = pd.read_excel ('Refixation_B.Nigra_HL_.xlsx') 
BnLL = pd.read_excel ('Refixation_B.Nigra_LL_.xlsx')
BrHL = pd.read_excel ('Refixation_B.Rapa_HL_.xlsx') 
BrLL = pd.read_excel ('Refixation_B.Rapa_LL_.xlsx')
AtHL = pd.read_excel ('Refixation_A.Thaliana_HL_.xlsx') 
AtLL = pd.read_excel ('Refixation_A.Thaliana_LL_.xlsx')

df=pd.DataFrame([])
df = df.append(HiHL); df = df.append(HiLL);
df = df.append(BnHL);df = df.append(BnLL) 
df = df.append(BrHL);df = df.append(BrLL);
df = df.append(AtHL);df = df.append(AtLL);

df.to_excel(PATH+'Refixation_all.xlsx')

Iinc=[75.02051333,100.0019467,124.9989333,149.9964667,175.0040333,
199.9953333,299.9946129,399.9886774,549.9749667,800.011,
1099.995333,1499.997097,1800.005333,2200.004194]

import seaborn as sns


fig, (ax1, ax2) = plt.subplots(1, 2)

ax1=sns.scatterplot(x=Iinc, y="f_refix",data=HiHL,label='Hi HL',s=48)
ax1=sns.scatterplot(x=Iinc, y="f_refix",data=BnHL,label='Bn HL',s=48)
ax1=sns.scatterplot(x=Iinc, y="f_refix",data=BrHL,label='Br HL',s=48)
ax1=sns.scatterplot(x=Iinc, y="f_refix",data=AtHL,label='At HL',s=48)
ax1.set_ylim(top=0.45)
ax1.set_ylim([0.3,0.45])
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

# plt.xlabel("Irradiance (µmol m${^2}$ s$^{-1}$)", fontsize= 24)
# plt.ylabel("f$_{refix}$ (-)", fontsize= 24)
# plt.ylabel("Fraction of refixed CO2 (-)", fontsize= 24)

ax2=sns.scatterplot(x=Iinc, y="f_refix",data=HiLL,label='Hi LL',s=48)
ax2=sns.scatterplot(x=Iinc, y="f_refix",data=BnLL,label='Bn LL',s=48)
ax2=sns.scatterplot(x=Iinc, y="f_refix",data=BrLL,label='Br LL',s=48)
ax2=sns.scatterplot(x=Iinc, y="f_refix",data=AtLL,label='At LL',s=48)
ax2.set_ylim([0.3,0.45])

plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("Irradiance (µmol m${^2}$ s$^{-1}$)", fontsize= 24)
plt.ylabel("Fraction of refixed CO${_2}$ (-)", fontsize= 24)

# plt.xlabel("Irradiance (µmol m${^2}$ s$^{-1}$)", fontsize= 24)
# plt.ylabel("f$_{refix}$ (-)", fontsize= 24)
# plt.ylabel("Fraction of refixed CO${_2}$ (-)", fontsize= 24)

plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

