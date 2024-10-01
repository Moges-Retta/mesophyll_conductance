# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 17:00:24 2021

@author: retta001
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
from tabulate import tabulate
import pandas as pd

import warnings
warnings.filterwarnings("ignore")


PATH = (r'\\WURNET.NL\Homes\retta001\My Documents\Project\2022\1-D Model\Resistance_to_CO2_diffusion_C3_leaf\\')


import tkinter as tk
from PIL import Image, ImageTk

root=tk.Tk()

canvas = tk.Canvas(root,width=600,height=300)
canvas.grid(columnspan=6,rowspan=6)
#ACI and AI
logo = Image.open('logo.tif')
logo = ImageTk.PhotoImage(logo)
logo_label = tk.Label(image=logo)
logo_label.image = logo


logo2 = Image.open('WUR_logo.png')
logo2 = ImageTk.PhotoImage(logo2)
logo_2_label = tk.Label(image=logo2)
logo_2_label.image = logo
logo_2_label.grid(column=4,row=3)

r = tk.IntVar()  # Mitochonria location, inner = 1, outer =0, everywhere 0<r<1
r.set("0")
r2 = tk.IntVar() # Gap effect, no change = 1, increase >1, decrease 0<=r2<1
r2.set("1")

r3 = tk.IntVar() # plants, H. incana = 1, B. nigra = 2
r3.set("1")

r4 = tk.IntVar() # treatment, HL = 1, LL=2
r4.set("1")

frame_mitochonria = tk.LabelFrame(root,text ="Mitochondria position relative to Chloroplast", padx=5,pady=5)
frame_mitochonria.grid(column=0,row=3)
frame_gap = tk.LabelFrame(root,text ="Chloroplast gap effect", padx=5,pady=5)
frame_gap.grid(column=0,row=2)
frame_species = tk.LabelFrame(root,text ="Plant", padx=5,pady=5)
frame_species.grid(column=0,row=0)
frame_treatments = tk.LabelFrame(root,text ="Light, High (HL) or low (LL) ", padx=5,pady=5)
frame_treatments.grid(column=0,row=1)

#mitochondria location and gap effect
btn1=tk.Radiobutton(frame_mitochonria,text="Inner",variable=r,value=1)
btn2=tk.Radiobutton(frame_mitochonria,text="Outer",variable=r,value=0)
btn3=tk.Radiobutton(frame_mitochonria,text="Everywhere",variable=r,value=3)
btn4=tk.Radiobutton(frame_gap,text="Increases",variable=r2,value=2)
btn5=tk.Radiobutton(frame_gap,text="Decreases",variable=r2,value=0)
btn6=tk.Radiobutton(frame_gap,text="No change",variable=r2,value=1)

#Spicies and treatment buttons
btn7=tk.Radiobutton(frame_species,text="H. incana",variable=r3,value=1)
btn8=tk.Radiobutton(frame_species,text="B. nigra",variable=r3,value=2)
btn11=tk.Radiobutton(frame_species,text="B. rapa",variable=r3,value=3)
btn12=tk.Radiobutton(frame_species,text="A. thaliana",variable=r3,value=4)
btn9=tk.Radiobutton(frame_treatments,text="HL",variable=r4,value=1)
btn10=tk.Radiobutton(frame_treatments,text="LL",variable=r4,value=2)


btn7.grid(column=0,row=0)
btn8.grid(column=1,row=0)
btn11.grid(column=0,row=1)
btn12.grid(column=1,row=1)
btn9.grid(column=0,row=0)
btn10.grid(column=1,row=0)

btn1.grid(column=0,row=2)
btn2.grid(column=1,row=2)
btn3.grid(column=2,row=2)

btn4.grid(column=0,row=3)
btn5.grid(column=1,row=3)
btn6.grid(column=2,row=3)

   
#Instructions
text_box = tk.Text(root,height=10,width=50)
text_box.grid(column=4,row=5,rowspan=2)
text_box.insert(1.0,"Compute the response of:\n \n1. photosynthesis to CO2 and light\n \n2. the fraction of re-fixed (photo)respired CO2\n\
                \n3. the sensitivity of photosynthesis to leaf anatomical parameters\n\
                \nModel is based on\n \nYin and Struik 2020, PS. 144(1). pp.85-99\
                \nBerghuijs et al. 2015,PS, 238,297-311\n \nProject: Extremophile\nMoges Retta")
text_box.tag_configure("center",justify="center")

physical_constants = Physical_constants()
[T,P,R,H,DCO2_water] = physical_constants.get_phyical_constants()

gap_effect = Names.GAP_NOCHANGE
mitochondria_location=Names.INNER_CYTOSOL
chloroplast_mitochondria_arrangement = [mitochondria_location.value,gap_effect.value]

treatment = 'HL'
plant = Plant(Names.HIncana.value,treatment,chloroplast_mitochondria_arrangement)

# Photosynthesis and refixation responses
O = 0.21
# gas_exch_measurement = Gas_exchange_measurement(O,incana.get_name(),treatment)    
#age = 15
gas_exch_measurement = Gas_exchange_measurement(O,plant,treatment)  

def get_plant(): 
    plant_code =  r3.get()
    treatment_code =  r4.get()
    location = r.get()
    gap = r2.get()
    if location==1:
        mitochondria_location=Names.INNER_CYTOSOL
    elif location==0:
        mitochondria_location=Names.OUTER_CYTSOL
    else:
        mitochondria_location=Names.EVERYWHERE
    
    if gap==1:
        gap_effect=Names.GAP_NOCHANGE
    elif gap>1:
        gap_effect=Names.GAP_INCREASES
    else:
        gap_effect=Names.GAP_DECREASES
    
    if treatment_code==1:
        treatment=Names.HL.value
    else:
        treatment=Names.LL.value
        
    if plant_code==1:
        plant_name=Names.HIncana.value
    elif plant_code==2:
        plant_name=Names.BNigra.value
    elif plant_code==3:
        plant_name=Names.BRapa.value
    elif plant_code==4:
        plant_name=Names.AThaliana.value
    
    
    chloroplast_mitochondria_arrangement = [mitochondria_location.value,gap_effect.value]
    
    plant = Plant(plant_name,treatment,chloroplast_mitochondria_arrangement)
    
    return plant

def sensitivity_object():
    plant = get_plant()
    O = 0.21
    gas_exch_measurement = Gas_exchange_measurement(O,plant.get_name(),plant.get_treatment())    
    sensitivity = Sensitivity(plant,gas_exch_measurement)
    return sensitivity

def calculate_response(curve,plant,I_ave,Ci_ave,O,gs_ave):
        f_refix_ci =[]; F_CI=[];f_refix_cell=[];f_refix_ias=[];
        AN_CI_mod=[];
                
        I =  [Iinc for Iinc in I_ave] 
        CI =  [Ci for Ci in Ci_ave] 
        df = pd.DataFrame([],columns=['Plant','Treatment','A','F','f_refix','f_refix_cell','f_refix_ias'])
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
        
        df['A']=np.asarray(AN_CI_mod).flatten()
        df['F']=np.asarray(F_CI).flatten()
        df['f_refix']=np.asarray(f_refix_ci).flatten()
        df['f_refix_cell']=np.asarray(f_refix_cell).flatten()
        df['f_refix_ias']=np.asarray(f_refix_ias).flatten()
        df.loc[:,'Plant']=plant.get_name()
        df.loc[:,'Treatment']=plant.get_treatment()
        df.to_excel(PATH + 'Refixation_'+plant.get_name()+'_'+plant.get_treatment()+ '_'+'.xlsx', index = False)
        return [AN_CI_mod,F_CI,f_refix_ci,f_refix_cell,f_refix_ias]
   
    
def plot_response(curve,x,y,y_mod,a_err):
        fig, ax = plt.subplots(figsize=(10, 10))
        plt.rcParams.update({'font.size': 16})
        plt.ylabel("Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)",fontsize=24)
        if curve==Names.CO2:
            ax.plot(x, y_mod,'k', label='Mod.',
                linewidth=1.4, markersize=10)
            ax.errorbar(x, y, a_err,fmt='ko',label='Expt.', markersize=10)
            plt.xlabel("Intercellular CO$_2$ (µmol mol$^{-1}$)",fontsize=24)  
        else:
            ax.plot(x, y_mod,'k', label='Mod.',linewidth=1.4, markersize=10)
            ax.errorbar(x, y, a_err,fmt='ko',label='Expt.', markersize=10)
            plt.xlabel("Irradiance (µmol m$^{-2}$ s$^{-1}$)",fontsize=24)  
            ax.tick_params(labelsize='medium', width=2)
            # ax.legend(loc='lower right', fontsize='x-large')
            
      
def calculate_values():
    plant = get_plant()
    plant2 = get_plant()
    O = 0.21
    gas_exch_measurement = Gas_exchange_measurement(O,plant.get_name(),plant.get_treatment())  

    # import copy
    # plant2 = copy.deepcopy(plant) # 
    
    parameters = Data(plant)
    
    df_r = pd.DataFrame([],columns=['Replicate','Tissue', 'Sm/S', 'Sc/s','rCO2','w'])
    df_r=[]
    df_r_i = pd.DataFrame([],columns=['Cell type','Cell component','Resistance','% tot'])
    replicate = [1,2,3]
    for r in replicate:
        
        [fi,ti_pal,peff_i,zeta_i,e_i] = parameters.get_palisade_wall(r)

        cell_wall = Cell_component(Names.CELLWALL,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i)
        myMesophyll = Make_mesophyll(parameters,cell_wall,"")
        [palisade,spongy] = myMesophyll.make_cells(r)
    
        """ calculate rdiff, An and w  """
        plant.add_cells(palisade)
        
        [rdiff_pali,*others] = plant.calculate_rliq()
        rdiff_pali=rdiff_pali*H/10**5
        Sc_S_pali = plant.calculate_Sc_S()
        Sm_S_pali = plant.calculate_Sm_S()
        
        plant.add_cells(spongy)
        plant2.add_cells(spongy)

        [rdiff_spongy,*others] = plant2.calculate_rliq()
        rdiff_spongy= rdiff_spongy*H/10**5
    
        [rliq,w] = plant.calculate_rliq()
        rliq= rliq*H/10**5
        rliq= rliq*(R*T/H) #equivalnet gas phase,R*T/H=0.83,

        rias = plant.calculate_rias(r)
        rias= rias*H/10**5  # H = Pa m3/mol, 10**5 Pa/bar == bar m2 s/mol
        
        rdiff = rias+ rliq
        
        Sc_S_tissue = plant.calculate_Sc_S()
        Sm_S_tissue = plant.calculate_Sm_S()
        
        Sc_S_spongy = Sc_S_tissue - Sc_S_pali
        Sm_S_spongy = Sm_S_tissue - Sm_S_pali
        
        df = plant.calculate_r_individual(r)
        df_r_i = df_r_i.append(df)
        df_pali = df[df['Cell type']=='PALISADE']
        w_pali = df_pali['w'].values[0]
        rdiff_pali = df_pali['Resistance'].values
        rdiff_pali = sum(rdiff_pali)
        df_spon = df[df['Cell type']=='SPONGY']
        w_spongy = df_spon['w'].values[0]   
        rdiff_spongy = df_spon['Resistance'].values
        rdiff_spongy = sum(rdiff_spongy)        
        
        r1 = [r,'Palisade',round(Sm_S_pali[0],2),round(Sc_S_pali[0],2),round(rdiff_pali,3),round(w_pali,3)]
        r2 = [r,'Spongy',round(Sm_S_spongy[0],2),round(Sc_S_spongy[0],2),round(rdiff_spongy,3),round(w_spongy,3)]
        r3 = [r,'IAS',0,0,round(rias,3),0]
        r4 = [r,'Liquid',0,0,round(rliq[0],3),0]
        r5 = [r,'Mesophyll',round(Sm_S_tissue[0],2),round(Sc_S_tissue[0],2),round(rdiff[0],3),round(w[0],3)]
        
        df_r.append(r1)
        df_r.append(r2)
        df_r.append(r3)
        df_r.append(r4)
        df_r.append(r5)
        
        
        plant = get_plant()
        plant2 = get_plant()
        

    df_r = pd.DataFrame(df_r,columns=['Replicate','Tissue', 'Sm/S', 'Sc/s','rCO2','w'])
    
    df_temp = df_r[df_r['Tissue']=='Palisade']
    df_temp=df_temp[['Sm/S', 'Sc/s','rCO2','w']]
    df_pali = df_temp.mean(axis=0)
    
    df_temp = df_r[df_r['Tissue']=='Spongy']
    df_temp=df_temp[['Sm/S', 'Sc/s','rCO2','w']]
    df_spo = df_temp.mean(axis=0)

    df_temp = df_r[df_r['Tissue']=='IAS']
    df_temp=df_temp[['Sm/S', 'Sc/s','rCO2','w']]
    df_ias = df_temp.mean(axis=0)
    
    df_temp = df_r[df_r['Tissue']=='Liquid']
    df_temp=df_temp[['Sm/S', 'Sc/s','rCO2','w']]
    df_liq = df_temp.mean(axis=0)
    
    df_temp = df_r[df_r['Tissue']=='Mesophyll']
    df_temp=df_temp[['Sm/S', 'Sc/s','rCO2','w']]
    df_mes = df_temp.mean(axis=0)
    
    Sm_S_pali = df_pali['Sm/S']
    Sc_S_pali = df_pali['Sc/s']
    rdiff_pali = df_pali['rCO2']
    w_pali = df_pali['w']
    
    Sm_S_spongy = df_spo['Sm/S']
    Sc_S_spongy = df_spo['Sc/s']
    rdiff_spongy = df_spo['rCO2']
    w_spongy = df_spo['w']
    
    Sm_S_tissue = df_mes['Sm/S']
    Sc_S_tissue = df_mes['Sc/s']
    
    rdiff = df_mes['rCO2']
    w = df_mes['w']
   
    rias = df_ias['rCO2']
    rliq = df_liq['rCO2']
    
    table = [['Tissue', 'Sm/S', 'Sc/s','rCO2','w'],
          ['Palisade', str(round(Sm_S_pali,2)), str(round(Sc_S_pali,2)),str(round(rdiff_pali,3)),str(round(w_pali,3))], 
          ['Spongy', str(round(Sm_S_spongy,2)),  str(round(Sc_S_spongy,2)),str(round(rdiff_spongy,3)),str(round(w_spongy,3))], 
          ['IAS',' ',' ',str(round(rias,3)),''], 
          ['Liquid',' ',' ',str(round(rliq,3)),''], 
          ['Mesophyll', str(round(Sm_S_tissue,2)), str(round(Sc_S_tissue,2)),str(round(rdiff,3)),str(round(w,3))]];
    print(tabulate(table, headers='firstrow', tablefmt='grid'))
    # df_r_i.to_excel(PATH + 'rCO2_individual_'+plant.get_name()+'_'+plant.get_treatment()+ '_'+'.xlsx', index = False)
    # df_r.to_excel(PATH + 'gm_anatomy_'+plant.get_name()+'_'+plant.get_treatment()+ '_'+'.xlsx', index = False)
    
    # df_r_i.to_excel(PATH + 'rCO2_individual_fit'+plant.get_name()+'_'+plant.get_treatment()+ '_'+'.xlsx', index = False)
    # df_r.to_excel(PATH + 'gm_anatomy_fit'+plant.get_name()+'_'+plant.get_treatment()+ '_'+'.xlsx', index = False)
    
    df_temp = df_r_i[df_r_i['Cell type']=='PALISADE']
    df_temp_cw = df_temp[df_temp['Cell component']=='CELLWALL']
    df_temp_pm = df_temp[df_temp['Cell component']=='PLASMALLEMA']
    df_temp_cy = df_temp[df_temp['Cell component']=='CYTOSOL']
    df_temp_ce = df_temp[df_temp['Cell component']=='CHLENVELOPE']
    df_temp_st = df_temp[df_temp['Cell component']=='STROMA']

    # df_temp_cw=df_temp_cw[['Resistance', '% tot']]
    df_temp_cw = df_temp_cw.mean(axis=0)

    # df_temp_ce=df_temp_ce[['Resistance', '% tot']]
    df_temp_ce = df_temp_ce.mean(axis=0)

    # df_temp_pm=df_temp_pm[['Resistance', '% tot']]
    df_temp_pm = df_temp_pm.mean(axis=0)
    
    # df_temp_cy=df_temp_cy[['Resistance', '% tot']]
    df_temp_cy = df_temp_cy.mean(axis=0)
    
    df_temp_st = df_temp_st.mean(axis=0)


    df=pd.DataFrame([])
    df=df.append(df_temp_cw,ignore_index=True)
    df=df.append(df_temp_pm,ignore_index=True)
    df=df.append(df_temp_ce,ignore_index=True)
    df=df.append(df_temp_cy,ignore_index=True)
    df=df.append(df_temp_st,ignore_index=True)

    df_temp = df_r_i[df_r_i['Cell type']=='SPONGY']
    df_temp_cw = df_temp[df_temp['Cell component']=='CELLWALL']
    df_temp_pm = df_temp[df_temp['Cell component']=='PLASMALLEMA']
    df_temp_cy = df_temp[df_temp['Cell component']=='CYTOSOL']
    df_temp_ce = df_temp[df_temp['Cell component']=='CHLENVELOPE']
    df_temp_st = df_temp[df_temp['Cell component']=='STROMA']


    # df_temp_cw=df_temp_cw[['Resistance', '% tot']]
    df_temp_cw = df_temp_cw.mean(axis=0)
    
    # df_temp_ce=df_temp_ce[['Resistance', '% tot']]
    df_temp_ce = df_temp_ce.mean(axis=0)
    
    # df_temp_pm=df_temp_pm[['Resistance', '% tot']]
    df_temp_pm = df_temp_pm.mean(axis=0)
    
    # df_temp_cy=df_temp_cy[['Resistance', '% tot']]
    df_temp_cy = df_temp_cy.mean(axis=0)
    
    df_temp_st = df_temp_st.mean(axis=0)

    df=df.append(df_temp_cw,ignore_index=True)
    df=df.append(df_temp_pm,ignore_index=True)
    df=df.append(df_temp_cy,ignore_index=True)
    df=df.append(df_temp_ce,ignore_index=True)
    df=df.append(df_temp_st,ignore_index=True)

    rtot = df['Resistance'].sum()
    df["% tot"]=round(df['Resistance']*100/rtot,2)
    
    df.loc[:,'Cell type']=['PALISADE','PALISADE','PALISADE','PALISADE','PALISADE','SPONGY','SPONGY','SPONGY','SPONGY','SPONGY']
    df.loc[:,'Cell component']=['CELLWALL','PLASMALLEMA','CYTOSOL','CHLENVELOPE','STROMA','CELLWALL','PLASMALLEMA','CYTOSOL','CHLENVELOPE','STROMA']
 
    df = pd.DataFrame(df,columns=['Cell type','Cell component','Resistance','% tot'])
    df=df.sort_values(by='% tot', ascending=False)
    print(tabulate(df, headers='keys', tablefmt='psql',showindex=False))  
    

    plant = get_plant()
    replicate = 4 # average leaf anatomical properties
    [fi,ti_pal,peff_i,zeta_i,e_i] = parameters.get_palisade_wall(replicate)

    cell_wall = Cell_component(Names.CELLWALL,fi,ti_pal,peff_i,zeta_i,DCO2_water,e_i)
    myMesophyll = Make_mesophyll(parameters,cell_wall,"")
    [palisade,spongy] = myMesophyll.make_cells(replicate)
    
    plant.add_cells(palisade)
    plant.add_cells(spongy)

    df_ave = gas_exch_measurement.average_A_CI()
    I_ave_ci = df_ave['Irradiance'].values
    Ci_ave_ci = df_ave['Intercellular_CO2_concentration'].values
    A_ave_ci = df_ave['Net_CO2_assimilation_rate'].values
    gs_ave_ci = df_ave['Stomatal_conductance_for_CO2'].values
    a_err =df_ave['Photo_err'].values   
                 
    [AN_CI_mod,F_CI,f_refix_ci,f_refix_cell_ci,f_refix_ias_ci] = calculate_response(Names.CO2, plant,I_ave_ci,Ci_ave_ci,O,gs_ave_ci)
                            
    plot_response(Names.CO2,Ci_ave_ci,A_ave_ci,np.asarray(AN_CI_mod).flatten(),a_err)
                            
    df_ave= gas_exch_measurement.average_A_I()
    I_ave_i = df_ave['Irradiance'].values
    Ci_ave_i = df_ave['Intercellular_CO2_concentration'].values
    A_ave_i = df_ave['Net_CO2_assimilation_rate'].values
    gs_ave_i = df_ave['Stomatal_conductance_for_CO2'].values
    a_err_i =df_ave['Photo_err'].values    
                       
    [AN_I_mod,F_I,f_refix_i,f_refix_cell_i,f_refix_ias_i] = calculate_response(Names.LIGHT,plant,I_ave_i,Ci_ave_i,O,gs_ave_i)
    plot_response(Names.LIGHT,I_ave_i,A_ave_i,np.asarray(AN_I_mod).flatten(),a_err_i)
        
    # Refixation
                            
    plt.rcParams["figure.figsize"] = (12,6)
    fig, ax = plt.subplots(1,2, sharey=True)
    ax[0].plot(Ci_ave_ci.transpose(),np.asarray(f_refix_ci)[0],'ko', label='Total',markersize=10)
    ax[0].plot((Ci_ave_ci.transpose()),np.asarray(f_refix_cell_ci)[0],'k^', label='Cell',fillstyle='none',markersize=10)
    ax[0].plot((Ci_ave_ci.transpose()),np.asarray(f_refix_ias_ci)[0],'ks', label='IAS',fillstyle='none',markersize=10)
    ax[0].set_ylabel("Fraction of refixation (-)",fontsize=24)
    ax[0].set_xlabel("Intercellular CO${_2}$ (µmol mol$^{-1}$)",fontsize=24)  
    ax[0].set_ylim(top=0.5)
    ax[1].plot( (I_ave_i),np.asarray(f_refix_i)[0],'ko', label='Total',markersize=10)
    ax[1].plot( (I_ave_i),np.asarray(f_refix_cell_i)[0],'k^', label='Cell',fillstyle='none',markersize=10)
    ax[1].plot( (I_ave_i),np.asarray(f_refix_ias_i)[0],'ks', label='IAS',fillstyle='none',markersize=10)
    ax[1].set_xlabel("Irradiance (µmol m$^{-2}$ s$^{-1}$)",fontsize=24)  
    ax[0].legend(loc='upper right',fontsize=14)
    ax[1].set_ylim(top=0.5)
                                    
    plt.tight_layout()
    plt.show()
    plt.rcParams["figure.figsize"] = (10,10)
        

""" Sensitivity """
def sensitivity_wall():
    sensitivity=sensitivity_object()
    cell_type_code =  r_s.get()
    sensitivity.of_cell_wall(cell_type_code)
    
    
def sensitivity_chloroplast_envelop():
    sensitivity=sensitivity_object()    
    cell_type_code =  r_s.get()  
    sensitivity.of_chloroplast_envelop(cell_type_code)
    
    
def sensitivity_cytosol():
    sensitivity=sensitivity_object()    
    cell_type_code =  r_s.get()      
    sensitivity.of_cytosol(cell_type_code)


def sensitivity_stroma():
    sensitivity=sensitivity_object()    
    cell_type_code =  r_s.get()      
    sensitivity.of_stroma(cell_type_code)
    

def sensitivity_palisade(): 
    sensitivity=sensitivity_object()    
    cell_type_code =  r_s.get() 
    sensitivity.of_palisade(cell_type_code)
    

def sensitivity_spongy(): 
    sensitivity=sensitivity_object()    
    cell_type_code =  r_s.get() 
    sensitivity.of_spongy(cell_type_code)


frame_compute = tk.LabelFrame(root,text ="Calculate photosynthesis")
frame_compute.grid(column=3,row=2)

frame_sensitivity = tk.LabelFrame(root,text ="Calculate sensitivity")
frame_sensitivity.grid(column=0,row=5)

r_s = tk.StringVar() # spo, pal
r_s.set(Names.PALISADE.value)

btn_s=tk.Radiobutton(frame_sensitivity,text="spongy",variable=r_s,value=Names.SPONGY.value)
btn_p=tk.Radiobutton(frame_sensitivity,text="palisade",variable=r_s,value=Names.PALISADE.value)

btn_s.grid(column=0,row=0)
btn_p.grid(column=1,row=0)


#Calculate rdiff
browse_text = tk.StringVar()
browse_btn = tk.Button(frame_compute,textvariable=browse_text, command=lambda:calculate_values(),font="Raleway",padx=25,pady=25,bd=5)
browse_text.set("COMPUTE")
browse_btn.grid(column=0,row=0)

#Calculate sensitivity
browse_2_text = tk.StringVar()
browse_2_btn = tk.Button(frame_sensitivity,textvariable=browse_2_text, command=lambda:sensitivity_wall(),font="Raleway",padx=25,pady=25,bd=5)
browse_2_text.set("Cell wall")
browse_2_btn.grid(column=0,row=2)

browse_3_text = tk.StringVar()
browse_3_btn = tk.Button(frame_sensitivity,textvariable=browse_3_text, command=lambda:sensitivity_chloroplast_envelop(),font="Raleway",padx=17,pady=25,bd=5)
browse_3_text.set("Chl.envlope")
browse_3_btn.grid(column=1,row=2)

browse_4_text = tk.StringVar()
browse_4_btn = tk.Button(frame_sensitivity,textvariable=browse_4_text, command=lambda:sensitivity_cytosol(),font="Raleway",padx=25,pady=25,bd=5)
browse_4_text.set("Cytosol")
browse_4_btn.grid(column=2,row=2)

browse_5_text = tk.StringVar()
browse_5_btn = tk.Button(frame_sensitivity,textvariable=browse_5_text, command=lambda:sensitivity_stroma(),font="Raleway",padx=27,pady=25,bd=5)
browse_5_text.set("Stroma")
browse_5_btn.grid(column=0,row=4)

browse_6_text = tk.StringVar()
browse_6_btn = tk.Button(frame_sensitivity,textvariable=browse_6_text, command=lambda:sensitivity_palisade(),font="Raleway",padx=26,pady=25,bd=5)
browse_6_text.set("Palisade")
browse_6_btn.grid(column=1,row=4)

browse_7_text = tk.StringVar()
browse_7_btn = tk.Button(frame_sensitivity,textvariable=browse_7_text, command=lambda:sensitivity_spongy(),font="Raleway",padx=25,pady=25,bd=5)
browse_7_text.set("Spongy")
browse_7_btn.grid(column=2,row=4)

root.after(2000)
root.mainloop()
