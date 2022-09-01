
import streamlit as st
import numpy as np
%matplotlib inline
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz as ctz

import plotly.graph_objects as go
from plotly.subplots import make_subplots

#----------------------------------------------

#concrete properties
Lf=0.1
Ccon=800
Kcon=1.7
ro=2200
Ccon=800

Lw=0.05
Af=10
Cw=Ccon*ro*Af*Lw
Rwcon=Lw/(Kcon*Af)

#----------------------------------------------

#air properties
Cp_air=1000
ro_air=1.2
V=30
C_air=ro_air*Cp_air*V

#slab film coefficient

hf=5
R12=1/(Af*hf)
Cf=Ccon*ro*Af*Lf
Cf1=Cf/2
Cf2=Cf1

#insulation for walls
Rins=2; #U<0.57
uw=1/Rins
r0w=2*1/(uw*Af)

#insulation for roof and floor

uff=1/4; #U<0.3
r0f=1/uff*Af;

#thermal resistance of the slab

Rf1f2=Lf/(2*Kcon*Af)
R2f1=Lf/(4*Kcon*Af)
Rf2g=R2f1

#h thermal resistance between air and side/back wall
R1sw1=1/(7*10)

#Window thermal resistance

Awin=6
uwin=1.5  #U<1.7
rwin=1/(uwin*Awin)
R1win=1/(5*Awin)

#Thermal resistance between indoor air and outdoor air including infiltration, windows and doors

R10=0.1
Delta_t=50

#------------------------------------------------------

T1=np.zeros((1200,1))
Tceiling=np.zeros((1200,1))
T2=np.zeros((1200,1))
Tf1=np.zeros((1200,1))
Tf2=np.zeros((1200,1))
Tsw1=np.zeros((1200,1))
Tsw2=np.zeros((1200,1))
Tbw=np.zeros((1200,1))
Twin=np.zeros((1200,1))

qaux=np.zeros((1200,1))
Error=np.zeros((1200,1))
PI=np.zeros((1200,1))

T1[0]=25

Tceiling[0]=25

Tf1[0]=25

Tf2[0]=25

Tsw1[0]=25

Tsw2[0]=25

Tbw[0]=25

Twin[0]=25


qaux[0]=0

qmax=2000

Kp=1000

Tg=20

import pandas as pd

Data= pd.read_csv('Data.csv')

Too_july=Data["To"]

QS_bw_july=Data["Qsolar_bw "]

QS_sw1_july=Data["Qsolar_ sw1"]

QS_sw2_july=Data["Qsolar_sw2"]

QS_window_july=Data["Qsolar_window"]

QS_ceil_july=Data["Qsolar_roof"]

Tsol_bw=Too_july+(0.3*QS_bw_july/10);

Tsol_sw1=Too_july+(0.3*QS_sw1_july/10);

Tsol_sw2=Too_july+(0.3*QS_sw2_july/10);

Tsol_ceil=Too_july+(0.3*QS_ceil_july/10);


Tsp=np.ones((1200,1))*24

Tsp2=np.ones((1200,1))*24

Tsp2[791:830]=20

Tsp2[1021:1060]=20

#----------------------------------------------------

pp=np.arange(0,1200-1)

for i in pp:
    
    if Tsp[i]<T1[i]:
        
        Error[i+1]=T1[i]-Tsp[i]
    else:
        Error[i+1]=0
    
    PI[i+1]=Kp*Error[i]
    if Error[i+1]>0.2:
        qaux[i+1]=PI[i+1]
    if PI[i+1]<qmax:
        qaux[i+1]=PI[i+1]
    if PI[i+1]>qmax:
            qaux[i+1]=qmax
   
    Delta_t=50
    
    T1[i+1]=(Delta_t/(15*C_air))*(((T2[i]-T1[i])/R12)+((Too_july[i]-T1[i])/R10)+((Tsw1[i]-T1[i])/R1sw1)+((Tsw2[i]-T1[i])/R1sw1)+((Tbw[i]-T1[i])/R1sw1)+((Twin[i]-T1[i])/R1win)+((Tceiling[i]-T1[i])/R1sw1)-qaux[i])+T1[i]
    
    # T1[i+1]=(Delta_t/(15*C_air))*(((T2[i]-T1[i])/R12)+((Too_july[i]-T1[i]/R10))+((Tsw1[i]-T1[i])/R1sw1)+((Tsw2[i]-T1[i])/R1sw1)+((Tbw[i]-T1[i])/R1sw1)+((Twin[i]-T1[i])/R1win)+((Tceiling[i]-T1[i])/R1sw1)-qaux[i])+T1[i] 
    
    Tsw1[i+1]=(Delta_t/Cw)*(((T1[i]-Tsw1[i])/R1sw1)+((Tsol_sw1[i]-Tsw1[i])/(r0w+Rwcon)))+Tsw1[i]
    
    Tsw2[i+1]=(Delta_t/Cw)*(((T1[i]-Tsw2[i])/R1sw1)+((Tsol_sw2[i]-Tsw2[i])/(r0w+Rwcon)))+Tsw2[i]
    
    Tbw[i+1]=(Delta_t/Cw)*(((T1[i]-Tbw[i])/R1sw1)+((Tsol_bw[i]-Tbw[i])/(r0w+Rwcon))+(0.3*0.5*QS_window_july[i]*10))+Tbw[i]
    
    Tceiling[i+1]=(Delta_t/Cw)*(((T1[i]-Tceiling[i])/R1sw1)+((Tsol_ceil[i]-Tceiling[i])/(r0w+r0f)))+Tceiling[i]

    Twin[i+1]=((T1[i]/R1sw1)+(Too_july[i]/rwin)+(0.1*QS_window_july[i]*6))/((1/R1sw1)+(1/rwin))
                                   
    T2[i+1]=((Tf1[i]/R2f1)+(T1[i]/R12)+(0.3*0.5*QS_window_july[i]*10))/((1/R2f1)+(1/R12))
                                   
    Tf1[i+1]=(Delta_t/Cf1)*(((Tf2[i]-Tf1[i])/Rf1f2)+((T2[i]-Tf1[i])/R2f1))+Tf1[i]
                                   
    Tf2[i+1]=(Delta_t/Cf2)*(((Tg-Tf2[i])/(Rf2g+r0f))+((Tf1[i]-Tf2[i])/Rf1f2))+Tf2[i]

    
    st.pyplot(plt.plot(qaux))

