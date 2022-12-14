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
st.set_page_config(layout="wide")

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
     value=(time(3, 00), time(8, 00)))

setpoint_bool = st.sidebar.checkbox('Flexible event for prev. days')

st.sidebar.markdown("""---""")

peak_time_select = st.sidebar.slider(
     "Peak time",
     value=(time(12, 00), time(18, 00))
     )



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
df_dr["COP"] = df_dr.apply(lambda x: chillerCOP(x["To"], x["qaux"]/qmax*100), axis = 1)
df_bs["COP"] = df_bs.apply(lambda x: chillerCOP(x["To"], x["qaux"]/qmax*100), axis = 1)

#Estimates electricity 
df_dr["e_w"] = df_dr["qaux"]/df_dr["COP"]
df_bs["e_w"] = df_bs["qaux"]/df_bs["COP"]


#Estimates KPIs
df_bs = df_bs.loc[date_day_str]
df_dr = df_dr.loc[date_day_str]
dx = 6*60

peak_hour_s = peak_time_select[0]
peak_hour_e = peak_time_select[1]

#Cooling KPIs--------------------------------------------------
    # Total energy used during the day
cooling_total_bs = np.trapz(df_bs['qaux'], dx=dx)/(1000*3600)
cooling_total_dr = np.trapz(df_dr['qaux'], dx=dx)/(1000*3600)

    # ADR event
hour_day_end = time(23,59)

df_sliced = df_bs.loc[hour_s:hour_e]
cooling_drevent_bs = np.trapz(df_sliced["qaux"], dx=dx)/(1000*3600)
df_sliced = df_dr.loc[hour_s:hour_e]
cooling_drevent_dr = np.trapz(df_sliced["qaux"], dx=dx)/(1000*3600)

cool_down_flex = (cooling_drevent_dr - cooling_drevent_bs)

    # After ADR event
df_sliced = df_bs.loc[hour_e:hour_day_end]
cooling_drafter_bs = np.trapz(df_sliced["qaux"], dx=dx)/(1000*3600)
df_sliced = df_dr.loc[hour_e:hour_day_end]
cooling_drafter_dr = np.trapz(df_sliced["qaux"], dx=dx)/(1000*3600)

cool_down_flex_after = (cooling_drafter_dr - cooling_drafter_bs)

    # Efficiency
cool_eff = cool_down_flex_after/cool_down_flex

    #Peak time analysis
bs_avgpeak_coolpower = df_bs.loc[peak_hour_s:peak_hour_e]["qaux"].mean()
dr_avgpeak_coolpower = df_dr.loc[peak_hour_s:peak_hour_e]["qaux"].mean()

cool_befi = (bs_avgpeak_coolpower - dr_avgpeak_coolpower)/1000
cool_befi_p = cool_befi/(bs_avgpeak_coolpower/1000)

    #Summary table
df_cool_table = pd.DataFrame({"CADR (kWh)":[cool_down_flex],
"Energy shift (kWh)":[cool_down_flex_after],
"E.Shift/CADR (%)":[cool_eff],
"Peak-h reduction (kW)": [cool_befi],
"Peak-h reduction (%)": [cool_befi_p]
})

    # Formats summary table
df_cool_table["CADR (kWh)"] = df_cool_table["CADR (kWh)"].map('{:,.2f}'.format)
df_cool_table["Energy shift (kWh)"] = df_cool_table["Energy shift (kWh)"].map('{:,.2f}'.format)
df_cool_table["E.Shift/CADR (%)"] = df_cool_table["E.Shift/CADR (%)"].map('{:,.2%}'.format)
df_cool_table["Peak-h reduction (kW)"] = df_cool_table["Peak-h reduction (kW)"].map('{:,.3f}'.format)
df_cool_table["Peak-h reduction (%)"] = df_cool_table["Peak-h reduction (%)"].map('{:,.1%}'.format)

#Electrical KPIs-----------------------------------------------
    # Total energy used during the day
elec_total_bs = np.trapz(df_bs['e_w'], dx=dx)/(1000*3600)
elec_total_dr = np.trapz(df_dr['e_w'], dx=dx)/(1000*3600)

    # ADR event
hour_day_end = time(23,59)

df_sliced = df_bs.loc[hour_s:hour_e]
elec_drevent_bs = np.trapz(df_sliced["e_w"], dx=dx)/(1000*3600)
df_sliced = df_dr.loc[hour_s:hour_e]
elec_drevent_dr = np.trapz(df_sliced["e_w"], dx=dx)/(1000*3600)

elec_down_flex = (elec_drevent_dr - elec_drevent_bs)

    # After ADR event
df_sliced = df_bs.loc[hour_e:hour_day_end]
elec_drafter_bs = np.trapz(df_sliced["e_w"], dx=dx)/(1000*3600)
df_sliced = df_dr.loc[hour_e:hour_day_end]
elec_drafter_dr = np.trapz(df_sliced["e_w"], dx=dx)/(1000*3600)

elec_down_flex_after = (elec_drafter_dr - elec_drafter_bs)

    # Efficiency
elec_eff = elec_down_flex_after/elec_down_flex

    #Peak time analysis
bs_avgpeak_elecpower = df_bs.loc[peak_hour_s:peak_hour_e]["e_w"].mean()
dr_avgpeak_elecpower = df_dr.loc[peak_hour_s:peak_hour_e]["e_w"].mean()

elec_befi = (bs_avgpeak_elecpower - dr_avgpeak_elecpower)/1000
elec_befi_p = elec_befi/(bs_avgpeak_elecpower/1000)

    # Summary table
df_elec_table = pd.DataFrame({"CADR (kWh)":[elec_down_flex],
"Energy shift (kWh)":[elec_down_flex_after],
"E.Shift/CADR (%)":[elec_eff],
"Peak-h reduction (kW)": [elec_befi],
"Peak-h reduction (%)": [elec_befi_p]})

    # Formats summary table
df_elec_table["CADR (kWh)"] = df_elec_table["CADR (kWh)"].map('{:,.2f}'.format)
df_elec_table["Energy shift (kWh)"] = df_elec_table["Energy shift (kWh)"].map('{:,.2f}'.format)
df_elec_table["E.Shift/CADR (%)"] = df_elec_table["E.Shift/CADR (%)"].map('{:,.2%}'.format)
df_elec_table["Peak-h reduction (kW)"] = df_elec_table["Peak-h reduction (kW)"].map('{:,.3f}'.format)
df_elec_table["Peak-h reduction (%)"] = df_elec_table["Peak-h reduction (%)"].map('{:,.1%}'.format)


#Plotting -------------------------------------------------------------------------------------   
fig_dr_day = make_subplots(specs=[[{"secondary_y": True}]])

df_bs = df_bs.loc[date_day_str]
df_dr = df_dr.loc[date_day_str]

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["e_w"],
    mode = 'lines',
    name = "HVAC Electricity baseline",
    line = dict(width = 1.0, color = "red", dash = "solid")
    ),secondary_y=False)

    #Adds metric
fig_dr_day.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["e_w"],
    mode = 'lines',
    name = "HVAC Electricity flexible",
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
    title="Electrical consumption",
    xaxis_title="Day Time",
    yaxis_title="HVAC electrical load (W)",
    legend_title="Variables",
    )

fig_dr_day.add_vline(x=datetime.combine(date_day_select,peak_hour_s),  line_width=1, line_dash="dot", line_color="grey")
fig_dr_day.add_vline(x=datetime.combine(date_day_select,peak_hour_e),  line_width=1, line_dash="dot", line_color="grey")

if sp != fstp:
        #Annotation CADR
    x1 = datetime.combine(date_day_select,hour_s).timestamp()
    x2 = datetime.combine(date_day_select,hour_e).timestamp()
    xa = (x2-x1)/2+x1
    CADR_x = datetime.fromtimestamp(xa)

    y2 = df_dr["e_w"].max()
    y1 = df_bs.loc[roundTime(CADR_x, 60*6)]["e_w"]
    ya = (y2-y1)*0.2+y1

    CADR_y = ya
    fig_dr_day.add_annotation(x=CADR_x, y=CADR_y,
                text="CADR = " + str(round(elec_down_flex,2)) + " kWh",
                showarrow=False,
                yshift=0)

        #Annotation Energy unloaded
    x1 = datetime.combine(date_day_select,hour_e).timestamp()
    x2 = datetime.combine(date_day_select,time(23,59)).timestamp()
    xe = (x2-x1)/2+x1
    energy_x = datetime.fromtimestamp(xe)

    y2 = df_bs["e_w"].max()
    y1 = df_dr["e_w"].min()
    ye = (y2-y1)*0.85+y1

    energy_y = ye
    fig_dr_day.add_annotation(x=energy_x, y=energy_y,
                text="Energy shift = " + str(round(elec_down_flex_after,2)) + " kWh",
                showarrow=False,
                yshift=0)


        #Annotation Ratio
    ratio_x = energy_x
    ratio_y = df_bs["e_w"].max()*1.25
    
    fig_dr_day.add_annotation(x=ratio_x, y=ratio_y,
                text="Energy shift/CADR ratio =" + str(round(elec_eff*100,2)) + "%",
                showarrow=False,
                yshift=0)

fig_dr_day.update_yaxes(title_text="Temperature (C)", secondary_y=True, titlefont = {"size": 20})
fig_dr_day.update_yaxes(secondary_y=False, titlefont = {"size": 20})
fig_dr_day.update_xaxes(titlefont = {"size": 20})


# Cooling graphs --------------------------------------------------------
fig_dr_day_cool = make_subplots(specs=[[{"secondary_y": True}]])

df_bs = df_bs.loc[date_day_str]
df_dr = df_dr.loc[date_day_str]

    #Adds metric
fig_dr_day_cool .add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["qaux"],
    mode = 'lines',
    name = "HVAC Cooling baseline",
    line = dict(width = 1.0, color = "orchid", dash = "solid")
    ),secondary_y=False)

    #Adds metric
fig_dr_day_cool .add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["qaux"],
    mode = 'lines',
    name = "HVAC Cooling flexible",
    line = dict(width = 2.0, color = "orchid", dash = "dash"),
    fill='tonexty'
    ),secondary_y=False)

    #Adds metric
fig_dr_day_cool .add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["tsp"],
    mode = 'lines',
    name = "Setpoint - baseline",
    line = dict(width = 1.0, color = "blue", dash='solid')
    ),secondary_y=True)

    #Adds metric
fig_dr_day_cool .add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["tsp2"],
    mode = 'lines',
    name = "Setpoint - flexible",
    line = dict(width = 2.0, color = "blue", dash='dash')
    ),secondary_y=True)

    #Adds metric
fig_dr_day_cool .add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 1.0, color = "green", dash='solid')
    ),secondary_y=True)

    #Adds metric
fig_dr_day_cool .add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["t1"],
    mode = 'lines',
    name = "Indoor air temperature",
    line = dict(width = 2.0, color = "green", dash='dash')
    ),secondary_y=True)


    #Adds metric
fig_dr_day_cool .add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["To"],
    mode = 'lines',
    name = "Ambient Temperature",
    line = dict(width = 1.0, color = "orange")
    ),secondary_y=True)

fig_dr_day_cool .update_layout(
    title="Cooling consumption",
    xaxis_title="Day Time",
    yaxis_title="HVAC cooling load (W)",
    legend_title="Variables",
    )

fig_dr_day_cool.add_vline(x=datetime.combine(date_day_select,peak_hour_s),  line_width=1, line_dash="dot", line_color="grey")
fig_dr_day_cool.add_vline(x=datetime.combine(date_day_select,peak_hour_e),  line_width=1, line_dash="dot", line_color="grey")

if sp != fstp:

        #Annotation CADR
    x1 = datetime.combine(date_day_select,hour_s).timestamp()
    x2 = datetime.combine(date_day_select,hour_e).timestamp()
    xa = (x2-x1)/2+x1
    CADR_x = datetime.fromtimestamp(xa)

    y2 = df_dr["qaux"].max()
    y1 = df_bs.loc[roundTime(CADR_x, 60*6)]["qaux"]
    ya = (y2-y1)*0.2+y1

    CADR_y = ya
    fig_dr_day_cool .add_annotation(x=CADR_x, y=CADR_y,
                text="CADR = " + str(round(cool_down_flex,2)) + " kWh",
                showarrow=False,
                yshift=0)

        #Annotation Energy unloaded
    x1 = datetime.combine(date_day_select,hour_e).timestamp()
    x2 = datetime.combine(date_day_select,time(23,59)).timestamp()
    xe = (x2-x1)/2+x1
    energy_x = datetime.fromtimestamp(xe)

    y2 = df_bs["qaux"].max()
    y1 = df_dr["qaux"].min()
    ye = (y2-y1)*0.85+y1

    energy_y = ye
    fig_dr_day_cool .add_annotation(x=energy_x, y=energy_y,
                text="Energy shift = " + str(round(cool_down_flex_after,2)) + " kWh",
                showarrow=False,
                yshift=0)


        #Annotation Ratio
    ratio_x = energy_x
    ratio_y = df_bs["qaux"].max()*1.25
    
    fig_dr_day_cool .add_annotation(x=ratio_x, y=ratio_y,
                text="Energy shift/CADR ratio =" + str(round(cool_eff*100,2)) + "%",
                showarrow=False,
                yshift=0)

fig_dr_day_cool .update_yaxes(title_text="Temperature (C)", secondary_y=True, titlefont = {"size": 20})
fig_dr_day_cool.update_yaxes(secondary_y=False, titlefont = {"size": 20})
fig_dr_day_cool.update_xaxes(titlefont = {"size": 20})


#Day selection plot-----------------------------------------------------------
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
fig_dr_year.add_vline(x=date_day_str,  line_width=1, line_dash="dash", line_color="red")
fig_dr_year.update_traces(marker_size=10)

fig_dr_year.update_layout(
    title="Weather on selected day",
    xaxis_title="Day Date",
    yaxis_title="Ambient Temperature (C)",
    legend_title="Variables",
    )

fig_dr_year.update_yaxes(secondary_y=False, titlefont = {"size": 20})
fig_dr_year.update_xaxes(titlefont = {"size": 20})


#Gauge plot cooling----------------------------------------------------------
bs_cons = cooling_total_bs
dr_cons = cooling_total_dr

if bs_cons >= dr_cons:
    color_gauge = "#3D9970"
else:
    color_gauge = "#FF4136"

fig_gauge_cool = go.Figure(go.Indicator(
    mode = "number+gauge+delta",
    value = dr_cons,
    number = {'valueformat':".2f"},
    domain = {'x': [0.1, 1], 'y': [0, 1]},
    title = {'text' :"<b>Cooling. (kWh)</b><br><span style='color: darkblue; font-size:0.7em'>Baseline " + str(round(bs_cons,2)) + " kWh</span>"},
    delta = {'reference': bs_cons, "relative" : True, 'valueformat':".2%", "increasing":{"color":"#FF4136", "symbol": "???"},"decreasing":{"color":"#3D9970", "symbol":"???"}},
    gauge = {
        'shape': "angular",
        'axis': {'range': [None, dr_cons*1.25]},
        'threshold': {
            'line': {'color': "darkblue", 'width': 3},
            'thickness': 1,
            'value': bs_cons},
        'steps': [
            {'range': [0, bs_cons], 'color': "orchid"},
            {'range': [bs_cons, dr_cons], 'color': color_gauge}],
        'bar':{'color':'indigo'}
        }))

fig_gauge_cool.update_layout(
    margin=dict(l=30, r=30, t=80, b=10),
    height = 250,
    width = 200
)

#Gauge plot electrical----------------------------------------------------------
bs_cons = elec_total_bs
dr_cons = elec_total_dr

if bs_cons >= dr_cons:
    color_gauge = "#3D9970"
else:
    color_gauge = "#FF4136"


fig_gauge_elec = go.Figure(go.Indicator(
    mode = "number+gauge+delta",
    value = dr_cons,
    number = {'valueformat':".2f"},
    domain = {'x': [0.1, 1], 'y': [0, 1]},
    title = {'text' :"<b>Elect. (kWh)</b><br><span style='color: darkblue; font-size:0.7em'>Baseline " + str(round(bs_cons,2)) + " kWh</span>"},
    delta = {'reference': bs_cons, "relative" : True, 'valueformat':".2%", "increasing":{"color":"#FF4136", "symbol": "???"},"decreasing":{"color":"#3D9970", "symbol":"???"}},
    gauge = {
        'shape': "angular",
        'axis': {'range': [None, dr_cons*1.25]},
        'threshold': {
            'line': {'color': "darkblue", 'width': 3},
            'thickness': 1,
            'value': bs_cons},
        'steps': [
            {'range': [0, bs_cons], 'color': "mistyrose"},
            {'range': [bs_cons, dr_cons], 'color': color_gauge}],
        'bar':{'color':'darkred'}
        }))


fig_gauge_elec.update_layout(
    margin=dict(l=30, r=30, t=80, b=10),
    height = 250,
    width = 200
)

# HVAC plots ---------------------------------------------------------------------------------
df_bs = df_bs.loc[date_day_str]
df_dr = df_dr.loc[date_day_str]

fig_CvsE = make_subplots(rows=1,
                         cols=2,
                         specs=[[{"secondary_y": True},{"secondary_y": True}]],
                         subplot_titles=("Loads vs COP", "COP vs Ambient temp")
                        )

#Adds metric
fig_CvsE.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["e_w"],
    mode = 'lines',
    name = "Electric load",
    line = dict(width = 2.0, color = "red", dash = "dash")
    ),row = 1, col =1,secondary_y=False)

#Adds metric
fig_CvsE.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["qaux"],
    mode = 'lines',
    name = "Cooling lod",
    line = dict(width = 2.0, color = "orchid", dash = "dash")
    ),row = 1, col =1,secondary_y=False)

#Adds metric
fig_CvsE.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["COP"],
    mode = 'lines',
    name = "COP (kW/kW)",
    line = dict(width = 3.0, color = "indigo", dash = "solid")
    ),row = 1, col =1,secondary_y=True)

#Adds metric
fig_CvsE.add_trace(go.Scatter(
    x=df_bs['date'],
    y=df_bs["COP"],
    mode = 'lines',
    name = "COP (kW/kW)",
    line = dict(width = 3.0, color = "indigo", dash = "solid")
    ),row = 1, col =2,secondary_y=False)

#Adds metric
fig_CvsE.add_trace(go.Scatter(
    x=df_dr['date'],
    y=df_dr["To"],
    mode = 'lines',
    name = "Ambient Temperature (C)",
    line = dict(width = 1.0, color = "orange", dash = "solid")
    ),row = 1, col =2,secondary_y=True)

fig_CvsE.update_layout(
    # title="HVAC efficiency",
    legend_title="Variables",
    )

fig_CvsE.update_yaxes(range=[0, 2000*1.1], title_text = "Power (W)", secondary_y = False, row =1, col =1)
fig_CvsE.update_yaxes(range=[0, 5], title_text = "COP (kW/kW)", secondary_y = True, row =1, col =1)
fig_CvsE.update_yaxes(range=[0, 5], title_text = "COP (kW/kW)", secondary_y = False, row =1, col =2)
fig_CvsE.update_yaxes(range=[0, 50], title_text = "Ambient Temperature (C)", secondary_y = True, row =1, col =2)

# Setting up page-----------------------------------------------------------------

    #To eras indices column
    # CSS to inject contained in a string
hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
                """


    #Title
st.markdown('<b style="color:darkgoldenrod ; font-size: 44px">Flexibility analysis for a simplified building in Dubai</b>', unsafe_allow_html=True)

st.markdown("""---""")
st.markdown('<b style="color:midnightblue ; font-size: 25px">Ambient and model conditions</b>', unsafe_allow_html=True)
col1, col2 = st.columns(2)
with col1:
   st.image("Box_model.jpg")

with col2:
   st.plotly_chart(fig_dr_year, use_container_width=True)

st.markdown("""---""")
st.markdown('<b style="color:midnightblue ; font-size: 25px">Model results</b>', unsafe_allow_html=True)

col1, col2 = st.columns([1,2])
with col1:
    st.plotly_chart(fig_gauge_cool, use_container_width=True)

    # Inject CSS with Markdown
    st.markdown(hide_table_row_index, unsafe_allow_html=True)
    st.table(df_cool_table)
    
with col2:
    st.plotly_chart(fig_dr_day_cool, use_container_width=True)


col3, col4 = st.columns([1,2])
with col3:
    st.plotly_chart(fig_gauge_elec, use_container_width=True)

    # Inject CSS with Markdown
    st.markdown(hide_table_row_index, unsafe_allow_html=True)
    st.table(df_elec_table)

with col4:
    st.plotly_chart(fig_dr_day, use_container_width=True)

st.markdown("""---""")
st.markdown('<b style="color:midnightblue ; font-size: 25px">HVAC efficiency</b>', unsafe_allow_html=True)
st.plotly_chart(fig_CvsE, use_container_width=True)











