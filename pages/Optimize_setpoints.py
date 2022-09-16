import streamlit as st

import pandas as pd
import numpy as np

from collections import defaultdict
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz as ctz
from function import simulate, download_data_csv, roundTime, chillerCOP

from datetime import datetime, timedelta, time

import plotly.graph_objects as go
from plotly.subplots import make_subplots


st.set_page_config(layout="wide")

sps = st.sidebar.slider(
     'Baseline setpoint',
     10, 35, (20,26))


# cost, df = sim_elec_cost_full(X_tsp, data_sim, date_day_str, cost_X)

#Sidebar
st.sidebar.image("gears.png")

    # Use inputs ---------------------------------------------------------------------------------
sps = st.sidebar.slider(
     'Baseline setpoint',
     10, 35, (20,26))

date_day_select = st.sidebar.date_input(
     "Analysis date",
     value = datetime(2019, 7, 1),
     min_value = datetime(2019, 1, 6),
     max_value = datetime(2019, 12, 30)
     )

setpoint_bool = st.sidebar.checkbox('Flexible event for prev. days')