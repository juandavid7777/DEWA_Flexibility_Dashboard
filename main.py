import pandas as pd
import numpy as np

from collections import defaultdict
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz as ctz
import streamlit as st
from function import simulate, download_data_csv, roundTime, chillerCOP

from datetime import datetime, timedelta, time

import plotly.graph_objects as go
from plotly.subplots import make_subplots

#Sidebar
st.sidebar.image("gears.png")

    # Use inputs ---------------------------------------------------------------------------------
sp = st.sidebar.slider(
     'Baseline setpoint',
     10, 35, 24)

fstp = st.sidebar.slider(
     'Flexible setpoint',
     10, 35, 24)

date_day_select = st.sidebar.date_input(
     "Analysis date",
     value = datetime(2019, 7, 1),
     min_value = datetime(2019, 1, 6),
     max_value = datetime(2019, 12, 30)
     )


dr_time = st.sidebar.slider(
     "Demand response event time:",
     value=(time(7, 00), time(13, 00)))

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
hour_s = dr_time[0]    #hour start
hour_e = dr_time[1]   #hour end

#Applies alternative setpoint
if setpoint_bool == False:
    data_sim.loc[date_day_str]["tsp2"] = data_sim.apply(lambda x: tsp2 if ((x["date"].time() >= hour_s) & (x["date"].time() < hour_e)) else tsp, axis = 1)

else:
    data_sim["tsp2"] = data_sim.apply(lambda x: tsp2 if ((x["date"].time() >= hour_s) & (x["date"].time() < hour_e)) else tsp, axis = 1)


# Simulation ----------------------------------------------------------------------------------
    #Gets parameters
param_names = ["Kp", "qmax", "Delta_t", "C_air", "Cw", "Cf1", "Cf2", "R12", "R10", "R1sw1", "R1win", "r0w", "Rwcon", "r0f", "rwin", "R2f1", "Rf1f2", "Rf2g", "Tg" ]
param_values = [Kp, qmax, Delta_t, C_air, Cw, Cf1, Cf2, R12, R10, R1sw1, R1win, r0w, Rwcon, r0f, rwin, R2f1, Rf1f2, Rf2g, Tg ]
params= dict(zip(param_names, param_values))

    #Set initialconditions
state_vars = ["t1", "tceil", "t2", "tf1", "tf2", "tsw1", "tsw2", "tbw", "twin", "qaux", "error", "pi"]
initial_values = [24,24,24,24,24,24,24,24,24,0,0,0]
x0 = dict(zip(state_vars, initial_values))

#Simulate baseline data
result = simulate(data_sim, params, x0, tsp = "tsp")
data_sim = data_sim.reset_index(drop = True)
df_bs = pd.concat([data_sim, result], axis=1).set_index("date", drop = False)

#Simulate flexible data
result = simulate(data_sim, params, x0, tsp = "tsp2")
data_sim = data_sim.reset_index(drop = True)
df_dr = pd.concat([data_sim, result], axis=1).set_index("date", drop = False)

#Estimates COP
df_dr["COP"] = df_dr.apply(lambda x: chillerCOP(x["To"], x["qaux"]/2000*100), axis = 1)
df_bs["COP"] = df_bs.apply(lambda x: chillerCOP(x["To"], x["qaux"]/2000*100), axis = 1)

#Estimates electricity 
df_dr["e_w"] = df_dr["qaux"]/df_dr["COP"]
df_bs["e_w"] = df_bs["qaux"]/df_bs["COP"]


#Estimates KPIs
    # Total energy used during the day
cooling_total_bs = np.trapz(df_bs['qaux'], dx=np.diff(df_bs['date'])/np.timedelta64(1, 's'))/(1000*3600)
cooling_total_dr = np.trapz(df_dr['qaux'], dx=np.diff(df_dr['date'])/np.timedelta64(1, 's'))/(1000*3600)

    # ADR event
hour_day_end = time(23,59)

df_sliced = df_bs.loc[hour_s:hour_e]
cooling_drevent_bs = np.trapz(df_sliced["qaux"], dx=np.diff(df_sliced['date'])/np.timedelta64(1, 's'))/(1000*3600)
df_sliced = df_dr.loc[hour_s:hour_e]
cooling_drevent_dr = np.trapz(df_sliced["qaux"], dx=np.diff(df_sliced['date'])/np.timedelta64(1, 's'))/(1000*3600)

down_flex = (cooling_drevent_dr - cooling_drevent_bs)

    # After ADR event
df_sliced = df_bs.loc[hour_e:hour_day_end]
cooling_drafter_bs = np.trapz(df_sliced["qaux"], dx=np.diff(df_sliced['date'])/np.timedelta64(1, 's'))/(1000*3600)
df_sliced = df_dr.loc[hour_e:hour_day_end]
cooling_drafter_dr = np.trapz(df_sliced["qaux"], dx=np.diff(df_sliced['date'])/np.timedelta64(1, 's'))/(1000*3600)

down_flex_after = (cooling_drafter_dr - cooling_drafter_bs)

    # Efficiency
eff = down_flex_after/down_flex



#Plotting -------------------------------------------------------------------------------------   
fig_dr_day = make_subplots(specs=[[{"secondary_y": True}]])

df_bs = df_bs.loc[date_day_str]
df_dr = df_dr.loc[date_day_str]

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["qaux"],
    mode = 'lines',
    name = "Cooling baseline",
    line = dict(width = 1.0, color = "red", dash = "solid")
    ),secondary_y=False)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["qaux"],
    mode = 'lines',
    name = "Cooling flexible",
    line = dict(width = 2.0, color = "red", dash = "dash"),
    fill='tonexty'
    ),secondary_y=False)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["tsp"],
    mode = 'lines',
    name = "Setpoint - baseline",
    line = dict(width = 1.0, color = "blue", dash='solid')
    ),secondary_y=True)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["tsp2"],
    mode = 'lines',
    name = "Setpoint - flexible",
    line = dict(width = 2.0, color = "blue", dash='dash')
    ),secondary_y=True)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 1.0, color = "green", dash='solid')
    ),secondary_y=True)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 2.0, color = "green", dash='dash')
    ),secondary_y=True)


    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["To"],
    mode = 'lines',
    name = "Ambient Temperature",
    line = dict(width = 1.0, color = "orange")
    ),secondary_y=True)

fig_dr_day.update_layout(
    title="Selection model results",
    xaxis_title="Time",
    yaxis_title="Cooling load (W)",
    legend_title="Variables",
    )

# if sp != fstp:

#         #Annotation CADR
#     x1 = datetime.combine(date_day_select,hour_s).timestamp()
#     x2 = datetime.combine(date_day_select,hour_e).timestamp()
#     xa = (x2-x1)/2+x1
#     CADR_x = datetime.fromtimestamp(xa)

#     y2 = df_dr["qaux"].max()
#     y1 = df_bs.loc[roundTime(CADR_x, 60*6)]["qaux"]
#     ya = (y2-y1)*0.2+y1

#     CADR_y = ya
#     fig_dr_day.add_annotation(x=CADR_x, y=CADR_y,
#                 text="CADR = " + str(round(down_flex,2)) + " kWh",
#                 showarrow=False,
#                 yshift=0)

#         #Annotation Energy unloaded
#     x1 = datetime.combine(date_day_select,hour_e).timestamp()
#     x2 = datetime.combine(date_day_select,time(23,59)).timestamp()
#     xe = (x2-x1)/2+x1
#     energy_x = datetime.fromtimestamp(xe)

#     y2 = df_bs["qaux"].max()
#     y1 = df_dr["qaux"].min()
#     ye = (y2-y1)*0.85+y1

#     energy_y = ye
#     fig_dr_day.add_annotation(x=energy_x, y=energy_y,
#                 text="Energy shift = " + str(round(down_flex_after,2)) + " kWh",
#                 showarrow=False,
#                 yshift=0)

#         #Annotation Ratio
#     ratio_x = energy_x
#     ratio_y = df_bs.loc[roundTime(ratio_x, 60*6)]["qaux"]*1.1
    
#     fig_dr_day.add_annotation(x=ratio_x, y=ratio_y,
#                 text="Energy shift/CADR ratio =" + str(round(eff*100,2)) + "%",
#                 showarrow=False,
#                 yshift=0)

fig_dr_day.update_yaxes(title_text="Temperature (C)", secondary_y=True)
fig_dr_day.update_yaxes(range=[0, 40], secondary_y = True)

#Second plot-----------------------------------------------------------
fig_dr_year = make_subplots(specs=[[{"secondary_y": True}]])

df_p = data

#Adds metric
fig_dr_year.add_trace(go.Scatter(
    x=df_p['date'],
    y=df_p["To"],
    mode = 'lines',
    name = "Ambient Temperature",
    line = dict(width = 1.0, color = "orange")
    ),secondary_y=False)

fig_dr_year.add_trace(go.Scatter(
    x=[date_day_str],
    y=[sp],
    mode = 'markers',
    name = "Setpoint selection",
    line = dict(width = 1.0, color = "red")
    ),secondary_y=False)

fig_dr_year.add_hline(y=sp,  line_width=1, line_dash="dash", line_color="blue")
#fig_dr_year.add_hline(y=fstp,  line_width=2, line_dash="dash", line_color="blue")
fig_dr_year.add_vline(x=date_day_str,  line_width=1, line_dash="dash", line_color="red")
fig_dr_year.update_traces(marker_size=10)

fig_dr_year.update_layout(
    title="Selection yearly overview",
    xaxis_title="Time",
    yaxis_title="Ambient Temperature (C)",
    legend_title="Variables",
    )


#Third plot----------------------------------------------------------

bs_cons = cooling_total_bs
dr_cons = cooling_total_dr


fig_gauge = go.Figure(go.Indicator(
    mode = "number+gauge+delta", value = dr_cons,
    domain = {'x': [0.1, 1], 'y': [0, 1]},
    title = {'text' :"<b>Cooling load</b><br><span style='color: darkgray; font-size:0.7em'>Baseline " + str(round(bs_cons,1)) + " kWh</span>"},
    delta = {'reference': bs_cons, "relative" : True},
    gauge = {
        'shape': "bullet",
        'axis': {'range': [None, dr_cons]},
        'threshold': {
            'line': {'color': "red", 'width': 2},
            'thickness': 0.75,
            'value': bs_cons},
        'steps': [
            {'range': [0, bs_cons], 'color': "lightgray"},
            {'range': [bs_cons, dr_cons], 'color': "lightsalmon"}]}))

fig_gauge.update_layout(height = 250)

# Setting up page

    #Title
st.markdown('<b style="color:darkgoldenrod ; font-size: 44px">Flexibility analysis for a simplified building in Dubai</b>', unsafe_allow_html=True)

    #Image
st.image("Box_model.jpg")
st.plotly_chart(fig_dr_year, use_container_width=True)
st.plotly_chart(fig_dr_day, use_container_width=True)
st.plotly_chart(fig_gauge, use_container_width=True)






