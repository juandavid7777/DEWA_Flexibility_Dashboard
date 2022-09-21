import streamlit as st

import pandas as pd
import numpy as np
import json

from collections import defaultdict
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz as ctz
from function import simulate, download_data_csv, roundTime, chillerCOP, sim_elec_cost_full

from datetime import datetime, timedelta, time

import plotly.graph_objects as go
from plotly.subplots import make_subplots


st.set_page_config(layout="wide")

#Sidebar
st.sidebar.image("gears.png")

    # Use inputs ---------------------------------------------------------------------------------
sps = st.sidebar.slider(
     'Comfort ranges',
     18, 26, (20,26),
     step = 2)

sp = sps[1]
fstp = sps[0]

date_day_select = st.sidebar.date_input(
     "Analysis date",
     value = datetime(2019, 7, 1),
     min_value = datetime(2019, 1, 6),
     max_value = datetime(2019, 12, 30)
     )

setpoint_bool = st.sidebar.checkbox('Flexible event for prev. days')

# Material properties -------------------------------------------------------------------------
Lf=0.1
Ccon=800
Kcon=1.7
ro=2200
Ccon=800

Lw=0.05
Af=10
Cw=Ccon*ro*Af*Lw
Rwcon=Lw/(Kcon*Af)

# Model properties -----------------------------------------------------------------------------
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
r0f=1/uff*Af

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

#Loads environmental data variables -------------------------------------------------------------
#Import time series data variables (U)
data = download_data_csv('environment_data.csv') 
data = data.set_index("date", drop = False)

#slices selected dates data for a n day simulation ---------------------------------------------
date_day_str = date_day_select.strftime("%Y-%m-%d")
date_day = date_day_select #datetime.strptime(date_day_str,'%Y-%m-%d')

d_5 = timedelta(days = 5)
d_p1 = timedelta(days = 2)

date_5dbefore = date_day-d_5
date_1dafter = date_day+d_p1

data_sim = data.loc[date_5dbefore:date_1dafter]

#Set point definitions
tsp = sp
data_sim["tsp"] = tsp
data_sim["tsp2"] = tsp

#Set point modification over period
tsp2 = fstp

#Cost profile - define a price profile
price_aed = 0.32*100
price_signal = download_data_csv('dynamic_setpoint_opt.csv') 
price_signal = price_signal.set_index("date", drop = False)
cost_X = list(price_signal["price_normalized"]*price_aed)


#Gets optimal setpoint

with open('optimized_dynamic_solutions.json', 'r') as fp:
    dict_sol = json.load(fp)
    
cost, df = sim_elec_cost_full(dict_sol[str(sp)+"-"+str(fstp)], data_sim, date_day_str, cost_X)

#Estimates COP
df["COP"] = df.apply(lambda x: chillerCOP(x["To"], x["qaux"]/2000*100), axis = 1)

#Estimates electricity 
df["e_w"] = df["qaux"]/df["COP"]

# First plot-------------------------------------------------------------------------------------   
fig_cost = make_subplots(specs=[[{"secondary_y": True}]])

    #Adds metric
fig_cost.add_trace(go.Scatter(
    x=df['date'],
    y=df["cost"],
    mode = 'lines',
    name = "Cost",
    line = dict(width = 1.0, color = "indigo", dash='solid')
    ),secondary_y=False)

    #Adds metric
fig_cost.add_trace(go.Scatter(
    x=df['date'],
    y=df["tsp2"],
    mode = 'lines',
    name = "Setpoint - flexible",
    line = dict(width = 2.0, color = "blue", dash='dash')
    ),secondary_y=True)

    #Adds metric
fig_cost.add_trace(go.Scatter(
    x=df['date'],
    y=df["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 2.0, color = "green", dash='dash')
    ),secondary_y=True)

#     #Adds metric
# fig_cost.add_trace(go.Scatter(
#     x=df['date'],
#     y=df["To"],
#     mode = 'lines',
#     name = "Ambient Temperature",
#     line = dict(width = 1.0, color = "orange")
#     ),secondary_y=True)

fig_cost.update_layout(
    # title="Model results",
    xaxis_title="Time",
    yaxis_title="HVAC electrical load (W)",
    legend_title="Variables",
    )

fig_cost.update_yaxes(title_text="Temperature (C)", secondary_y=True)
# fig_cost.update_yaxes(range=[0, 65], secondary_y = True)

# Second plot ---------------------------------------------------------

fig_dr_day = make_subplots(specs=[[{"secondary_y": True}]])

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df['date'],
    y=df["qaux"],
    mode = 'lines',
    name = "Cooling load",
    line = dict(width = 1.0, color = "red", dash = "solid")
    ),secondary_y=False)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df['date'],
    y=df["tsp2"],
    mode = 'lines',
    name = "Setpoint - flexible",
    line = dict(width = 2.0, color = "blue", dash='dash')
    ),secondary_y=True)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df['date'],
    y=df["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 2.0, color = "green", dash='dash')
    ),secondary_y=True)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df['date'],
    y=df["To"],
    mode = 'lines',
    name = "Ambient Temperature",
    line = dict(width = 1.0, color = "orange")
    ),secondary_y=True)

fig_dr_day.update_layout(
    # title="Model results",
    xaxis_title="Time",
    yaxis_title="HVAC electrical load (W)",
    legend_title="Variables",
    )

fig_dr_day.update_yaxes(title_text="Temperature (C)", secondary_y=True)
# fig_dr_day.update_yaxes(range=[0, 65], secondary_y = True)

# Layout --------------------------------------------------------------
st.plotly_chart(fig_cost, use_container_width=True)
st.plotly_chart(fig_dr_day, use_container_width=True)