import pandas as pd
import numpy as np

from collections import defaultdict
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz as ctz

from datetime import datetime,timedelta

import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Use inputs ---------------------------------

fstp = st.slider(
     'Select a flexible set point',
     15, 35, 18)
st.write('Values:', fstp)




# Material properties --------------------------
Lf=0.1
Ccon=800
Kcon=1.7
ro=2200
Ccon=800

Lw=0.05
Af=10
Cw=Ccon*ro*Af*Lw
Rwcon=Lw/(Kcon*Af)

# Model properties ----------------------------

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

#Insulation for walls
Rins=2; #U<0.57
uw=1/Rins
r0w=2*1/(uw*Af)

#Insulation for roof and floor
uff=1/4; #U<0.3
r0f=1/uff*Af;

#Thermal resistance of the slab
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

#Maximum capacity of the HVAC system
qmax=2000

#Proportional controler constant
Kp=1000

#Ground temperature
Tg=20

#Loads environmental data variables ---------------------------

#Import time series data variables (U)
data = pd.read_csv('Data.csv', parse_dates = ["date"], dayfirst  =True)
data = data.set_index("date", drop = False)

#Estimates solar temperatures
solar_vars = ['Qsolar_window', 'Qsolar_bw', 'Qsolar_roof', 'Qsolar_sw1', 'Qsolar_sw2',]

for var in solar_vars:
    
    suffix = var.split("_")
    var_name = "Tsol_" + suffix[1]
    data[var_name] = data["To"] + (0.3* data[var]/10)


#slices selected dates data for a n day simulation --------------
date_day_str = "2019-07-01"
date_day = datetime.strptime(date_day_str,'%Y-%m-%d')

d_5 = timedelta(days = 5)
d_p1 = timedelta(days = 2)

date_5dbefore = date_day-d_5
date_1dafter = date_day+d_p1

data_sim = data.loc[date_5dbefore:date_1dafter]

#Set point definitions
tsp =24
data_sim["tsp"] = tsp
data_sim["tsp2"] = tsp

#Set point modification over period
tsp2 = 20
hour_s = 7    #hour start
hour_e = 12   #hour end

#Applies alternative setpoint
data_sim.loc[date_day_str]["tsp2"] = data_sim.apply(lambda x: tsp2 if ((x["date"].hour >= hour_s) & (x["date"].hour < hour_e)) else tsp, axis = 1)
data_sim["tsp2"].plot() #.loc[date_day_str]

# Simulation
#Set initialconditions
state_vars = ["t1", "tceil", "t2", "tf1", "tf2", "tsw1", "tsw2", "tbw", "twin", "qaux", "error", "pi"]
initial_values = [24,24,24,24,24,24,24,24,24,0,0,0]
x0 = dict(zip(state_vars, initial_values))

def simulate(data,x0,tsp = "tsp2"):

    stats = defaultdict(list)

    #Runs simulation
    for i in range(1,len(data)+1):

        #if else for the i-1 values depending on initial conditions
        if i == 1:
            t1_1 = x0["t1"]
            tceil_1 = x0["tceil"]
            t2_1 = x0["t2"]
            tf1_1 = x0["tf1"]
            tf2_1 = x0["tf2"]
            tsw1_1 = x0["tsw1"]
            tsw2_1 = x0["tsw2"]
            tbw_1 = x0["tbw"]
            twin_1 = x0["twin"]
            qaux_1 = x0["qaux"]
            error_1 = x0["error"]
            pi_1 = x0["pi"]

        else:
            t1_1 = t1
            tceil_1 = tceil
            t2_1 = t2
            tf1_1 = tf1
            tf2_1 = tf2
            tsw1_1 = tsw1
            tsw2_1 = tsw2
            tbw_1 = tbw
            twin_1 = twin
            qaux_1 = qaux
            error_1 = error
            pi_1 = pi

        #Stores values for reporting
        stats['t1'].append(t1_1)
        stats['tceil'].append(tceil_1)
        stats['t2'].append(t2_1)
        stats['tf1'].append(tf1_1)
        stats['tf2'].append(tf2_1)
        stats['tsw1'].append(tsw1_1)
        stats['tsw2'].append(tsw2_1)
        stats['tbw'].append(tbw_1)
        stats['twin'].append(twin_1)
        stats['qaux'].append(qaux_1)
        stats['error'].append(error_1)
        stats['pi'].append(pi_1)
        
        #Controller modeling
            #PI estimation
        pi = Kp * error_1

            #Error calculation
        if data.iloc[i-1][tsp]<t1_1:
            error = t1_1-data[tsp][i-1]

        else:
            error = 0

            #Conditionals
        if error>0.2:
            qaux=pi

        if (0 < pi) & (pi<qmax):
            qaux=pi
        else:
            qaux = 0

        if pi>qmax:
            qaux=qmax

        #Estimation of next time step
        t1    = (Delta_t/(15*C_air))*(((t2_1-t1_1)/R12)+((data.iloc[i-1]["To"]-t1_1)/R10)+((tsw1_1-t1_1)/R1sw1)+((tsw2_1-t1_1)/R1sw1)+((tbw_1-t1_1)/R1sw1)+((twin_1-t1_1)/R1win)+((tceil_1-t1_1)/R1sw1)-qaux_1)+t1_1  
        tsw1  = (Delta_t/Cw)*(((t1_1-tsw1_1)/R1sw1)+((data.iloc[i-1]["Tsol_sw1"]-tsw1_1)/(r0w+Rwcon)))+tsw1_1
        tsw2  = (Delta_t/Cw)*(((t1_1-tsw2_1)/R1sw1)+((data.iloc[i-1]["Tsol_sw2"]-tsw2_1)/(r0w+Rwcon)))+tsw2_1
        tbw   = (Delta_t/Cw)*(((t1_1-tbw_1)/R1sw1)+((data.iloc[i-1]["Tsol_bw"]-tbw_1)/(r0w+Rwcon))+(0.3*0.5*data.iloc[i-1]["Qsolar_window"]*10))+tbw_1
        tceil = (Delta_t/Cw)*(((t1_1-tceil_1)/R1sw1)+((data.iloc[i-1]["Tsol_roof"]-tceil_1)/(r0w+r0f)))+tceil_1
        twin  = ((t1_1/R1sw1)+(data.iloc[i-1]["To"]/rwin)+(0.1*data.iloc[i-1]["Qsolar_window"]*6))/((1/R1sw1)+(1/rwin))                              
        t2    = ((tf1_1/R2f1)+(t1_1/R12)+(0.3*0.5*data.iloc[i-1]["Qsolar_window"]*10))/((1/R2f1)+(1/R12))                              
        tf1   = (Delta_t/Cf1)*(((tf2_1-tf1_1)/Rf1f2)+((t2_1-tf1_1)/R2f1))+tf1_1                              
        tf2   = (Delta_t/Cf2)*(((Tg-tf2_1)/(Rf2g+r0f))+((tf1_1-tf2_1)/Rf1f2))+tf2_1

    #Results into dataframe
    result = pd.DataFrame().from_dict(stats)
    
    return result

#Simulate baseline data
result = simulate(data_sim, x0, tsp = "tsp")
data_sim = data_sim.reset_index(drop = True)
df_bs = pd.concat([data_sim, result], axis=1).set_index("date", drop = False)


#Simulate baseline data
result = simulate(data_sim, x0, tsp = "tsp2")
data_sim = data_sim.reset_index(drop = True)
df_dr = pd.concat([data_sim, result], axis=1).set_index("date", drop = False)

    
fig = make_subplots(specs=[[{"secondary_y": True}]])

df_bs = df_bs.loc[date_day_str]
df_dr = df_dr.loc[date_day_str]

#Adds metric
fig.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["qaux"],
    mode = 'lines',
    name = "Cooling baseline",
    line = dict(width = 1.0, color = "red", dash = "solid")
    ),secondary_y=False)

#Adds metric
fig.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["qaux"],
    mode = 'lines',
    name = "Cooling flexible",
    line = dict(width = 2.0, color = "red", dash = "dash")
    ),secondary_y=False)

#Adds metric
fig.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["tsp"],
    mode = 'lines',
    name = "Setpoint - baseline",
    line = dict(width = 1.0, color = "blue", dash='solid')
    ),secondary_y=True)

#Adds metric
fig.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["tsp2"],
    mode = 'lines',
    name = "Setpoint - flexible",
    line = dict(width = 2.0, color = "blue", dash='dash')
    ),secondary_y=True)

#Adds metric
fig.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 1.0, color = "green", dash='solid')
    ),secondary_y=True)

#Adds metric
fig.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 2.0, color = "green", dash='dash')
    ),secondary_y=True)


#Adds metric
fig.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["To"],
    mode = 'lines',
    name = "Ambient Temperature",
    line = dict(width = 1.0, color = "orange")
    ),secondary_y=True)

fig.update_yaxes(range=[0, 40], secondary_y = True)

st.plotly_chart(fig)

