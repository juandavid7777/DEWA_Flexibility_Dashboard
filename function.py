from collections import defaultdict
import pandas as pd
from datetime import datetime, timedelta, time
import streamlit as st

@st.cache(suppress_st_warning=True)
def simulate(data,param, x0, tsp = "tsp"):

    stats = defaultdict(list)

    #Gets model parameters
    Kp = param["Kp"]
    qmax = param["qmax"]
    Delta_t = param["Delta_t"]
    C_air = param["C_air"]
    Cw = param["Cw"]
    Cf1 = param["Cf1"]
    Cf2 = param["Cf2"]
    R12 = param["R12"]
    R10 = param["R10"]
    R1sw1 = param["R1sw1"]
    R1win = param["R1win"]
    r0w = param["r0w"]
    Rwcon = param["Rwcon"]
    r0f = param["r0f"]
    rwin = param["rwin"]
    R2f1 = param["R2f1"]
    Rf1f2 = param["Rf1f2"]
    Rf2g = param["Rf2g"]
    Tg = param["Tg"]

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

@st.cache(suppress_st_warning=True)
def download_data_csv(url):

    df = pd.read_csv(url, parse_dates = ["date"], dayfirst  =True)

    return df


def roundTime(dt=None, roundTo=60):
    """ Round a datetime object to any time lapse in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
    if dt == None : 
        dt = datetime.datetime.now()
    seconds = (dt.replace(tzinfo=None) - dt.min).seconds
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + timedelta(0,rounding-seconds,-dt.microsecond)

@st.cache(suppress_st_warning=True)
def chillerCOP(to, loading):
    x = to
    y = loading
    
    C1 = -0.0000357229471371693
    C2 = 4.29222732259951E-06
    C3 = 6.36621617952334E-07
    C4 = 0.00001150137344614
    C5 = 0.00403204029230074
    C6 = -0.001466186389158
    C7 = -0.00126546932281316
    C8 = -0.225123563022733
    C9 = 0.114233895801329
    C10 = 6.66022445845032
    
    COP = (C1*x**3) + (C2*y**3) + (C3*x**2*y) + (C4*x*y**2) + (C5*x**2) + (C6*y**2) + (C7*x*y) + (C8*x) + (C9*y) + C10
    
    return COP


@st.cache(suppress_st_warning=True)
def sim_elec_cost_full(tsp_X, data_sim, date_day_str, cost_X):

    #Param values

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
    
    #Gets parameters
    param_names = ["Kp", "qmax", "Delta_t", "C_air", "Cw", "Cf1", "Cf2", "R12", "R10", "R1sw1", "R1win", "r0w", "Rwcon", "r0f", "rwin", "R2f1", "Rf1f2", "Rf2g", "Tg" ]
    param_values = [Kp, qmax, Delta_t, C_air, Cw, Cf1, Cf2, R12, R10, R1sw1, R1win, r0w, Rwcon, r0f, rwin, R2f1, Rf1f2, Rf2g, Tg ]
    params= dict(zip(param_names, param_values))

    #Set initial conditions
    state_vars = ["t1", "tceil", "t2", "tf1", "tf2", "tsw1", "tsw2", "tbw", "twin", "qaux", "error", "pi"]
    initial_values = [24,24,24,24,24,24,24,24,24,0,0,0]
    x0 = dict(zip(state_vars, initial_values))

    #Replaces setpoint
    data_sim = data_sim.set_index("date", drop = False)
    data_sim.loc[date_day_str, "tsp2"] = tsp_X

    #Simulate baseline data
    result = simulate(data_sim, params, x0, tsp = "tsp2")
    data_sim = data_sim.reset_index(drop = True)
    df_bs = pd.concat([data_sim, result], axis=1).set_index("date", drop = False)

    df_bs = df_bs.loc[date_day_str]
    df_bs["cost"] = cost_X
    df_bs["coste"] = df_bs["cost"] * df_bs['qaux']

    cost_total = df_bs["coste"].sum()
    
    return cost_total, df_bs    