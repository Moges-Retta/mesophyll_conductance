# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 12:47:15 2021

@author: angel
"""
import math
# cell wall
fi = 1.0;
ti_pal = 0.120*10**-6;
ti_spo = 0.117*10**-6;
peff_i = 0.20;
zeta_i = 1.0;
DCO2_water = 1.79*10**-9 ;                 
H = 2941;
gamma_tissue_pali = 1.497                
Lm_L_tissue_pali = 6.01    
f_pal = 0.5140184418014165
Lm_L_tissue_spongy = 6.27
gamma_tissue_spongy= 1.281    


ti=ti_pal*f_pal+(1-f_pal)*ti_spo

Ri = ti*fi/(peff_i*zeta_i*DCO2_water)

Sm_s = Lm_L_tissue_pali*gamma_tissue_pali + Lm_L_tissue_spongy*gamma_tissue_spongy

Rl = Ri*H/10**5/Sm_s

""" photosynthesis"""
Jmax = 263.7  
alpha2LL= 0.721        
teta = 0.760
k2LL = 0.381          
phi2 = 0.721 
s = 0.526

Vcmax = 256
KmC = 267             
KmO = 164        
Rd = 2.46
Sco = 3260
Tp = 12.6 

alpha = 0

H = 2941; 
Iinc_d= [1499.62,
1499.87,
1499.82,
1499.7,
1499.6,
1498.88,
1498.98,
1501.28,
1500.11,
1500.29,
1500.21,
]

Ci_d = [250.761,
195.782,
138.333,
84.9497,
56.9371,
431.692,
603.863,
783.236,
960.881,
1327.73,
1685.62]

O = 210
gama_star = 0.5*O/Sco
rdiff = 3.94
w = 0.67
    
def calculate_J(Iinc):
    A = k2LL*Iinc+Jmax;
    B = 4*teta*Jmax*k2LL*Iinc;
    return (A-math.sqrt(A**2-B))/(2*teta)

def limiting_A(case,Iinc,Ci):
    if case=="Aj":
        x1=calculate_J(Iinc)/4
        x2=2*gama_star
    elif case=="Ac":
        x1=Vcmax
        x2=KmC*(1+O/KmO)
    else:
        x1=3*Tp
        x2=-gama_star
    a = x2+gama_star*(1-w*(1-alpha))
    b = w*(1-alpha)*(Rd*x2+gama_star*x1)-(x2+gama_star*(1-w*(1-alpha)))*(x1-Rd)-(x2+gama_star)*(Ci+x2)/rdiff
    
            b = m*(Rd*x2+gama_star*x1)-(x2+gama_star*(1-m))*(x1-Rd)-(photosynthesis.get_Ci()+x2)*((x2+gama_star)/photosynthesis.get_rdiff()+sigma*(x1-Rd))
        
    c = -w*(1-alpha)*(Rd*x2+gama_star*x1)*(x1-Rd)+(x2+gama_star)*(x1*(Ci-gama_star)-Rd*(Ci+x2))/rdiff

    if b**2-4*a*c<0:
        A=3*Tp-Rd
    else:
        A = (-b-math.sqrt(b**2-4*a*c))/(2*a)
    return A

i = 0
AN_CI_mod=[]
for Iinc in Iinc_d:
    Ci = Ci_d[i];
    Ac = limiting_A("Ac",Iinc,Ci)
    Aj = limiting_A("Aj",Iinc,Ci)
    Ap = limiting_A("Ap",Iinc,Ci)
    A = min(Ac,Aj,Ap)
    AN_CI_mod.append(A)
    i+=1


""" Based on Yin&Struik 2020"""
from Electron import Electron
photosynthesis = Photosynthesis(admiro,Iinc,Ci)  

electron = Electron(Jmax,k2LL,teta,phi2,photosynthesis.get_Iinc())
J=electron.calculated_J()
        # putative transition point based on  Yin&Struik 2020,
        #10.1007/s11120-020-00716-z
Cc_ts = max((24*gama_star*Tp+(1+3*aa)*gama_star*J)/(J-12*Tp),(1.0001+3*aa)*gama_star)
m = photosynthesis.get_w()*(1-alpha)
t1 = (Cc_ts-(1-m)*gama_star)*3*Tp
t2 = Cc_ts-(1+3*aa)*gama_star
gm_ts = 1/photosynthesis.get_rdiff() + sigma*(3*Tp/t2);
Ci_ts = Cc_ts+(t1/t2-(1-m)*Rd)/gm_ts
#        print(self.Ci)
#        print(Ci_ts)
   if photosynthesis.get_Ci()>Ci_ts:
        x1=J/4
        x2=2*gama_star
        [a,b,c]= photosynthesis.coefficients(x1,x2)
        AJ = (-b-math.sqrt(b**2-4*a*c))/(2*a)
    #            print(AJ)
        x1=Vcmax
        x2=KmC*(1+O/KmO)
        [a,b,c]= photosynthesis.coefficients(x1,x2)
        AC = (-b-math.sqrt(b**2-4*a*c))/(2*a)
    #            print(-b+math.sqrt(b**2-4*a*c))
        x1=3*Tp
        x2=-(1+3*aa)*gama_star
        [a,b,c]= photosynthesis.coefficients(x1,x2)
        if b**2-4*a*c<0:
            AP = min(AC,AJ)
        else:
            AP = (-b-math.sqrt(b**2-4*a*c))/(2*a)
    A = min(AC,AJ,AP)
    else:
        x1=J/4
        x2=2*gama_star
        [a,b,c]= photosynthesis.coefficients(x1,x2)
        AJ = (-b-math.sqrt(b**2-4*a*c))/(2*a)
        x1=Vcmax
        x2=KmC*(1+O/KmO)
        [a,b,c]= photosynthesis.coefficients(x1,x2)
        AC = (-b-math.sqrt(b**2-4*a*c))/(2*a)
    A = min(AC,AJ)
    
        
    def coefficients(self,x1,x2):
        m = photosynthesis.get_w()*(1-alpha)
        a = x2+gama_star*(1-m)+sigma*(photosynthesis.get_Ci()+x2)
        b = m*(Rd*x2+gama_star*x1)-(x2+gama_star*(1-m))*(x1-Rd)
        -(photosynthesis.get_Ci()+x2)*((x2+gama_star)/photosynthesis.get_rdiff()+sigma*(x1-Rd))
        -sigma*((x1*(photosynthesis.get_Ci()-gama_star)-Rd*(photosynthesis.get_Ci()+x2)))
        c = -m*(Rd*x2+gama_star*x1)*(x1-Rd)
        +((x2+gama_star)/photosynthesis.get_rdiff()+sigma*(x1-Rd))*(x1*(photosynthesis.get_Ci()-gama_star)
        -Rd*(photosynthesis.get_Ci()+x2))
        return [a,b,c]
        