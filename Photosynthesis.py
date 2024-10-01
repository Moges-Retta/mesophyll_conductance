# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:33:32 2021

@author: Moges
                  
"""
import math
from Physical_constants import Physical_constants
import numpy as np
import pandas as pd
from Data import Data

physical_constants = Physical_constants()
[T,P,R,H,DCO2_water] = physical_constants.get_phyical_constants()  
# Model based on Yin and Struik 2018. 
#Model becomes that of Berghuijs et al when alpha,aa and sigma are,zero,         
aa = 0.3 # assume α = 0.3, the fraction of glycolate car- bon not returned to the chloroplast
sigma = 0 # sigma 1.41 for tomatao, phenomenological model of gm from Yin&Struik 2020

class Photosynthesis:
    """   
            Model Photosynthesis 
            An :  net photosynthesis (µmol m-2 s-1)
            
            An model is based on gm model presented in  Berghuijs et al. 2015, PS,10.1016/j.plantsci.2015.06.022 
            in additon to that of An Yin&Struik 10.1007/s11120-017-0340-8
            
            lambda = 1 (all mito- chondria stay closely behind chloroplasts as if carboxylation 
            and (photo)respiratory CO2 release occur in one compartment)
            
            lambda = 0, corresponding to the original model of Tholen et al. (2012) that 
            applies to the case where all mitochondria locate in the outer cytosol.
            
            Jmax            # Maximum rate of electron transport through Photosystem II 
                              at saturating light, 
                              µmol e− m−2 leaf s−1
                              
            teta            # Convexity factor of the response of J to Iinc                  
            
                              
            alpha2LL        # Quantum yield of electron transport through Photosystem II
                              under strictly electron-transport-limiting conditions on the
                              basis of light absorbed by both Photosystem I and Photosystem
                              II,
                              mol e− mol−1 photon
            
            phi2            # Quantum yield of electron transport through Photosystem II, 
                              mol e− mol−1 photon
                              
            k2LL            # Conversion factor of incident irradiance into electron 
                              transport under electron-transport-limited conditions, 
                              mol e− mol−1 photon
            
            s               # Slope of the assumed linear relationship between AN and 1
                            # 1/4*Iinc*Phi2 under strictly electron-transport-limited conditions
                            
            KmC             # Michaelis–Menten constant of Rubisco for CO2, 
                              µbar CO2
            
            KmO             # Michaelis–Menten constant of Rubisco for O, 
                              mbar
            Tp              # Rate of triose phosphate utilization, 
                              µmol phosphate m−2 leaf s−1
            
                            
            Sco             # Relative CO2/O2 specificity factor of Rubisco, 
                              mbar O2 µbar−1 CO2    
                              
            gamma_star      # CO2 compensation point, 
                              µbar CO2
            
            lambda_mitochondria
            
            rsc             The fraction of mitochondria in the inner cytosol
            
            k               Fqctor allowing an increase(k>1), no change (k=0) and decrease (0<=k<1)
                            in the fraction of inner photorespired CO2 caused by chloroplast gaps
            
            rcx             Carboxylation resistance
            
            rdiff           # Total resistance for CO2 transport of the physical barriers 
                              in the mesophyll, 
                              m2 leaf s bar CO2 mol−1 CO2
                              
            rch            # Lumped resistance for CO2 transport of the chloroplast 
                              envelope and the stroma, and half the resistance of the 
                              cytosol, 
                              m2 leaf s bar CO2 mol−1 CO2
                              
            rwp             # Lumped resistance for CO2 transport of the cell wall, 
                              the plasmamembrane, and half the resistance of the cytosol,
                              m2 leaf s bar CO2 mol−1 CO2
    """
    
    def __init__(self,plant,Iinc,Ci,O,gs):
        self.plant=plant
        self.set_Iinc(Iinc)
        self.set_Ci(Ci)
        self.set_O(O)
        self.gs = gs
        # self.rdiff=plant.calculate_rliq()[0]*H/10**5+plant.calculate_rias()*H/10**5*(R*T/H)/plant.calculate_Sm_S()
        # self.w=plant.calculate_rliq()[1]
        
       
           
            
    def set_Iinc(self,Iinc):
        if Iinc>=0:
            self.Iinc = Iinc
            
            
    def set_Ci(self,Ci):
        if Ci>=0:
            self.Ci = Ci*P/10**5
            
            
    # def set_rdiff(self,rdiff):
    #     if rdiff>=0:
    #         self.rdiff = rdiff  
            
            
    # def set_w(self,w):
    #     if w>=0:
    #         self.w = w    
            
            
    def set_O(self,O):
        if O == 210 or O == 20: # ubar
            self.O = O
            
            
        if O==0.21 or O==0.020:
            self.O = O*1000
    
    
    def get_Iinc(self):
        return self.Iinc
    
    
    def get_Ci(self):
        return self.Ci
    
    
    def get_O(self):
        return self.O
    
    
    # def get_rdiff(self):
    #     return self.rdiff
    
    
    # def get_w(self):
    #     return self.w
    
    
    def get_A_Ci_mod(self):
        return self.A_Ci_mod
    
    
    def get_A_I_mod(self):
        return self.A_I_mod
    
    
    def coefficients(self,x1,x2):
        [Vcmax,KmC,KmO,Rd,Sco,Tp] = self.plant.get_rubiscoConstants()
        gama_star = 0.5*self.get_O()/Sco
        alpha = self.plant.get_alpha()

        m = self.get_w()*(1-alpha)
        a = x2+gama_star*(1-m)+sigma*(self.get_Ci()+x2)
        b = m*(Rd*x2+gama_star*x1)-(x2+gama_star*(1-m))*(x1-Rd)-(self.get_Ci()+x2)*((x2+gama_star)/self.get_rdiff()+sigma*(x1-Rd))-sigma*((x1*(self.get_Ci()-gama_star)-Rd*(self.get_Ci()+x2)))
        c = -m*(Rd*x2+gama_star*x1)*(x1-Rd)+((x2+gama_star)/self.get_rdiff()+sigma*(x1-Rd))*(x1*(self.get_Ci()-gama_star)-Rd*(self.get_Ci()+x2))
        return [a,b,c]
    
   
    def calculated_J(self):
        thisPlant = self.plant
        A = thisPlant.get_k2LL()*self.get_Iinc()+thisPlant.get_Jmax()
        B = (4*thisPlant.get_teta()*thisPlant.get_Jmax()*thisPlant.get_k2LL()*self.get_Iinc())
        return (A-math.sqrt(A**2-B))/(2*thisPlant.get_teta())
    
    def calculate_Cc(self,m,x1,x2,Rd,sigma,gama_star,Ci,gm0):
        A = -m*(Rd*x2+gama_star*x1)
        B = gm0*(x2 +gama_star)+sigma*(x1-Rd)
        C = x1*(Ci-gama_star)-Rd*(Ci+x2)
        return A*(x1-Rd)+B*C
    
    
    def calculate_A_FvCB(self):
        """ Calculate photosynthesis rate.
        
        Parameters:
        -----------------
        
        arg_1 : Object
        Plant object
        
        Returns:
        ----------------
        Array of floats

                          
        """            
        
        thisPlant = self.plant
        [Vcmax,KmC,KmO,Rd,Sco,Tp] = thisPlant.get_rubiscoConstants()
        [Jmax,teta,k2LL]  = thisPlant.get_eTransportConstants()
        gama_star = 0.5*self.get_O()/Sco
        J= self.calculated_J()
        alpha = self.plant.get_alpha()
        [rdiff,w] = self.plant.calculate_rdiff()
        rdiff= rdiff*H/10**5
        # putative transition point based on  Yin&Struik 2020,
        #10.1007/s11120-020-00716-z
        Cc_ts = max((24*gama_star*Tp+(1+3*aa)*gama_star*J)/(J-12*Tp),(1.0001+3*aa)*gama_star)
        m = self.get_w()*(1-alpha)
        t1 = (Cc_ts-(1-m)*gama_star)*3*Tp
        t2 = Cc_ts-(1+3*aa)*gama_star
        gm_ts = 1/self.get_rdiff() + sigma*(3*Tp/t2);
        Ci_ts = Cc_ts+(t1/t2-(1-m)*Rd)/gm_ts
        Ci = self.get_Ci()
       # AJ
        x1=J/4
        x2=2*gama_star
        [a,b,c]= self.coefficients(x1,x2)
        AJ = (-b-math.sqrt(b**2-4*a*c))/(2*a)
        CcJ = (gama_star*x1+x2*(AJ+Rd))/(x1-(AJ+Rd))
        FJ = gama_star*x1/(CcJ+x2)
        rcx_J = (CcJ+x2)/x1
        
        # AC
        x1=Vcmax
        x2=KmC*(1+self.get_O()/KmO)
        [a,b,c]= self.coefficients(x1,x2)
        AC = (-b-math.sqrt(b**2-4*a*c))/(2*a)
        CcC = (gama_star*x1+x2*(AC+Rd))/(x1-(AC+Rd))
        FC = gama_star*x1/(CcC+x2)
        rcx_C = (CcC+x2)/x1
        
        # AP
        x1=3*Tp
        x2=-(1+3*aa)*gama_star
        [a,b,c]= self.coefficients(x1,x2)            
        AP = (-b+math.sqrt(max(0.00001,b**2-4*a*c)))/(2*a) 
        CcP = (gama_star*x1+x2*(AP+Rd))/(x1-(AP+Rd))
        FP = gama_star*x1/(CcP+x2)
        rcx_P = (CcP+x2)/x1
        test_T = b**2-4*a*c;
            
        if Ci<Ci_ts or test_T<0:
            A = min(AC,AJ)
            F = min(FJ,FC)
            Cc = min(CcC,CcJ)
            rcx = min(rcx_C,rcx_J)
        else:
            A = min(AC,AJ,AP)
            F = min(FJ,FC,FP)
            Cc = min(CcC,CcJ,CcP)
            rcx = min(rcx_C,rcx_J,rcx_P)
        return [A,F,Cc,rcx]
    
    
    def calculate_f_refix(self): ## Eq. S3.4 of Yin & Struik 2020
        [*others,rcx] = self.calculate_A()
        rsc = 1*P/10**5/(self.gs) # m2 s bar/mol
        alpha = self.plant.get_alpha()
        rdiff=self.plant.calculate_rliq()[0]*H/10**5
        w=self.plant.calculate_rliq()[1]

        r_refix = alpha/rcx + (1-alpha)/(w*rdiff + rcx)
        r_escape = alpha/rcx + (1-alpha)/(w*rdiff + rcx) + alpha/(rdiff + rsc)+(1-alpha)/((1-w)*rdiff + rsc)
        f_refix = r_refix/(r_refix+r_escape)
        r_escape = alpha/rcx + (1-alpha)/(w*rdiff + rcx) + alpha/rdiff + (1-alpha)/((1-w)*rdiff)  ## rsc = 0
        f_refix_cell = r_refix/(r_refix+r_escape)
        f_refix_ias = f_refix-f_refix_cell
        

        return [f_refix,f_refix_cell,f_refix_ias]
    
    def calculate_A(self):
        physical_constants = Physical_constants()
        [T,P,R_const,H,DCO2_water] = physical_constants.get_phyical_constants()  
        """
        FvCB model of photosynthesis extended to include nitrogen assimilation
        Model based on Yin et al. 2021 10.1111/pce.14070 
        
        """
        thisPlant = self.plant
        df_vcmax = thisPlant.get_rubiscoConstants()
        df_jmax = thisPlant.get_eTransportConstants()

        w=thisPlant.calculate_rliq()[1]

        alpha = self.plant.get_alpha()
        m = w*(1-alpha);      

        kmc = df_vcmax['KmC'].unique()            
        kmo = df_vcmax['KmO'].unique()         
        sco = df_vcmax['Sco'].unique()   

        ci = self.get_Ci()
        O = self.get_O()
        if O==0.21:
            O=210
        i_inc = self.get_Iinc()
        alphaT = 0;
        alphaG = 0;
        alphaS=df_vcmax['alphaS'].unique() 

        # R = self.plant.get_sigma();
        R = 0;

        
        replicates = df_vcmax['Replicate'].unique()
        df_r = pd.DataFrame([],columns=['Amod','Cc_ave','F_ave','rcx_ave'])
        for replicate in replicates:
            df_vcmax_r=df_vcmax[df_vcmax['Replicate']==replicate];
            df_jmax_r=df_jmax[df_jmax['Replicate']==replicate]
            theta = df_jmax_r['theta'].values.astype(float)  
            rd = df_vcmax_r['Rd'].values.astype(float)
            jmax = df_jmax_r['Jmax'].values.astype(float)
            vcmax = df_vcmax_r['Vcmax'].values.astype(float)
            TP = df_vcmax_r['Tp'].values.astype(float)
            k2LL = df_jmax_r['k2LL'].values.astype(float)  
            
            rliq = thisPlant.calculate_rliq()[0];
            rias = thisPlant.calculate_rias(replicate);
            rdiff = rias + (R_const*T/H)*rliq[0];
            
            rdiff = rdiff*H/10**5; # bar m2 s/mol
            
            if R==0:
                gm0 = 1/rdiff;                
            else:
                gm0 = 0;
                R = self.plant.get_sigma();
                
            # Rubisco-limited part

            gamma_st = 0.5*O/sco;
            gamma_x = gamma_st*(1-alphaG+2*alphaT);
            
            
            #  Rubisco-limited part;
            KMCMO = kmc*(1+O/kmo);
            X1R = vcmax;
            X2R = KMCMO;
            
            PR = gm0*(X2R+gamma_x)+(X1R-rd)*R;
            QR2 = (ci-gamma_x)*X1R-(ci+X2R)*rd;
            QR = (X2R+gamma_x*(1-m))*(X1R-rd);
            AAR = X2R+gamma_x*(1-m)+R*(ci+X2R);
            BBR = m*(rd*X2R+gamma_x*X1R)-QR-PR*(ci+X2R)-R*QR2;
            CCR = -m*(rd*X2R+gamma_x*X1R)*(X1R-rd)+PR*QR2;
            AR = (-BBR-(BBR**2-4.*AAR*CCR)**0.5)/(2.*AAR);
            CcC = (gamma_st*X1R+X2R*(AR+rd))/(X1R-(AR+rd))
            FC = gamma_st*X1R/(CcC+X2R)
            rcx_C = (CcC+X2R)/X1R

            # Electron transport limited part;
            
            BB = k2LL*i_inc + jmax;
            J = (BB-(BB**2-4*theta*jmax*k2LL*i_inc)**0.5)/(2*theta);
            X1R = J/4;
            X2R = (1+2*alphaG-alphaT+alphaS)*2.*gamma_st;
            PR = gm0*(X2R+gamma_x)+(X1R-rd)*R;
            QR2 = (ci-gamma_x)*X1R-(ci+X2R)*rd;
            QR = (X2R+gamma_x*(1-m))*(X1R-rd);
            AAR = X2R+gamma_x*(1-m)+R*(ci+X2R);
            BBR = m*(rd*X2R+gamma_x*X1R)-QR-PR*(ci+X2R)-R*QR2;
            CCR = -m*(rd*X2R+gamma_x*X1R)*(X1R-rd)+PR*QR2;
            AJ = (-BBR-(BBR**2-4.*AAR*CCR)**0.5)/(2.*AAR);
            
            CcJ = (gamma_st*X1R+X2R*(AJ+rd))/(X1R-(AJ+rd))
            FJ = gamma_st*X1R/(CcJ+X2R)
            rcx_J = (CcJ+X2R)/X1R
            
            #  TPU limited part;            
            X1R = 3*TP;
            X2R = -(1+3*alphaG+6*alphaT+4*alphaS)*gamma_st;
            PR = gm0*(X2R+gamma_x)+(X1R-rd)*R;
            QR2 = (ci-gamma_x)*X1R-(ci+X2R)*rd;
            QR = (X2R+gamma_x*(1-m))*(X1R-rd);
            AAR = X2R+gamma_x*(1-m)+R*(ci+X2R);
            BBR = m*(rd*X2R+gamma_x*X1R)-QR-PR*(ci+X2R)-R*QR2;
            CCR = -m*(rd*X2R+gamma_x*X1R)*(X1R-rd)+PR*QR2;

            #NOTE: + sign in front of SQRT part;
            bacT = BBR**2-4*AAR*CCR;
            if AAR==0:  # when alphaG,alphaT and alphas = 0
                AP = X1R;
            else:
                AP = (-BBR + (np.maximum(0.00001,BBR**2-4*AAR*CCR))**0.5)/(2*AAR);
            
            #--- Finding the Ci transition point from Aj to Ap;
            CCTS = np.maximum((1.0001+3*alphaG+6*alphaT+4*alphaS)*gamma_x, (24*gamma_x*TP+(1+3*alphaG+6*alphaT+4*alphaS)*gamma_x*J)/(J-12*TP));
            GMTS = gm0+R*3*TP/(CCTS-(1+3*alphaG+6*alphaT+4*alphaS)*gamma_x);
            AATS = (CCTS-(1-m)*gamma_x)*3*TP/(CCTS-(1+3*alphaG+6*alphaT+4*alphaS)*gamma_x);
            CITS  = CCTS + (AATS-(1-m)*rd)/GMTS;

            #-- Model for net photosynthesis rate;
            if ci<CITS or bacT<0:
                A_mod = np.minimum(AR,AJ)
                F = np.minimum(FJ,FC)
                Cc = np.minimum(CcC,CcJ)
                rcx = np.minimum(rcx_C,rcx_J)
            else:
                A1 = np.minimum(AR,AJ)
                A_mod = np.minimum(A1,AP)
                CcP = (gamma_st*X1R+X2R*(AP+rd))/(X1R-(AP+rd))
                FP = gamma_st*X1R/(CcP+X2R)
                rcx_P = (CcP+X2R)/X1R
                FJ = np.array(FP, dtype=float)
                FC = np.array(FP, dtype=float)
                FP = np.array(FP, dtype=float)

                CcC = np.array(FP, dtype=float)
                CcJ = np.array(FP, dtype=float)
                CcP = np.array(FP, dtype=float)
                
                rcx_C = np.array(rcx_C, dtype=float)
                rcx_J = np.array(rcx_J, dtype=float)
                rcx_P = np.array(rcx_P, dtype=float)
                
                F = np.minimum(FJ,FC,FP)
                Cc = np.minimum(CcC,CcJ,CcP)
                rcx = np.minimum(rcx_C,rcx_J,rcx_P)

            df_r.loc[replicate-1,'Amod']=A_mod 
            df_r.loc[replicate-1,'Cc_ave']=Cc
            df_r.loc[replicate-1,'F_ave']=F
            df_r.loc[replicate-1,'rcx_ave']=rcx
            
        A_ave= np.mean(df_r['Amod'].values)
                    
        Cc_ave= np.mean(df_r['Cc_ave'].values)
                    
        F_ave= np.mean(df_r['F_ave'].values)
                    
        rcx_ave= np.mean(df_r['rcx_ave'].values)
        df_r = pd.DataFrame([],columns=['Amod','Cc_ave','F_ave','rcx_ave'])
        return [A_ave,F_ave,Cc_ave,rcx_ave]
