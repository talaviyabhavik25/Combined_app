import streamlit as st
import pandas as pd
import numpy as np

@st.cache_data
def cache_slip_table():
    url = 'http://raw.githubusercontent.com/Ahmedhassan676/Combined_app/main/positive_pump/slip.csv'
    
    return pd.read_csv(url)
df_slip = cache_slip_table()
def main_positive():
    html_temp="""
    <div style="background-color:lightblue;padding:16px">
    <h2 style="color:black"; text-align:center> Positive pump capacity calculations </h2>
    </div>
    <style>
    table {
    font-family: arial, sans-serif;
    border-collapse: collapse;
    width: 100%;
    }

    td, th {
    border: 1px solid #dddddd;
    text-align: left;
    padding: 8px;
    }


    </style>


        """
    st.markdown(html_temp, unsafe_allow_html=True)
    s1 = st.selectbox('Single or double acting?',('Single Acting','Double Acting'), key = 'type')
    
    d = st.number_input('Insert piston/plunger diameter (mm)', key = 'd')
    rpm = st.number_input('Insert motor speed (rpm)', key = 'rpm')
    num = st.selectbox('Insert  number of pistons/plungers',(1,2,3,4), key = 'number_pistons')
    st_len = st.number_input('Insert stroke length (mm)', key = 'st_len')
    if s1 == 'Single Acting':
        if st.button("Reveal pump capacity"):
                try:
                    A = np.pi*np.power(d,2)*0.25
                    Q = A*rpm*num*st_len*0.00006
                    st.success('Your pump capacity is {} m3/hr ({} L/hr) assuming no slip'.format(Q,Q*1000))
                    s2 = st.selectbox('correct pump slip based on viscosity?',('No','Yes'), key = 'slip')
                    if s2 == 'Yes':
                         vis = st.number_input('Insert fluid viscosity (C.st)', key = 'vis')
                         slip=np.interp(vis,df_slip['Viscosity (Centistokes)'],df_slip['Slip']) 
                         Q_corrected = 0.01*Q*slip
                         st.success('Your pump capacity is {} m3/hr ({} L/hr) assuming no slip'.format(Q_corrected,Q_corrected*1000))
                except (ValueError, TypeError, UnboundLocalError): st.write('Please Check your dataset(s)')

    else: 
        a = st.number_input('Insert rod diameter (mm)', key = 'rod_diameter') 
        try:
                    A = np.pi*np.power(d,2)*0.25
                    Q = (2*A-a)*rpm*num*st_len*0.00006
                    st.success('Your pump capacity is {} m3/hr ({} L/hr) assuming no slip'.format(Q,Q*1000))
                    s2 = st.selectbox('correct pump slip based on viscosity?',('No','Yes'), key = 'slip')
                    if s2 == 'Yes':
                         vis = st.number_input('Insert fluid viscosity (C.st)', key = 'vis')
                         slip=np.interp(vis,df_slip['Viscosity (Centistokes)'],df_slip['Slip']) 
                         Q_corrected = 0.01*Q*slip
                         st.success('Your pump capacity is {} m3/hr ({} L/hr) assuming no slip'.format(Q_corrected,Q_corrected*1000))
        except (ValueError, TypeError, UnboundLocalError): st.write('Please Check your dataset(s)')


if __name__ == '__main__':
    main_positive()