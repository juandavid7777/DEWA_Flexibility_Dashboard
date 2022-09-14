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

@st.cache
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

@st.cache    
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