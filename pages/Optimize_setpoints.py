import streamlit as st

st.write("Hello Juan")

sp = st.sidebar.slider(
     'Baseline setpoint',
     10, 35, (20,26))