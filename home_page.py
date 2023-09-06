import streamlit_extras as st_ext
from streamlit_extras import colored_header
import streamlit as st
from streamlit_extras.switch_page_button import switch_page

description = 'This tool was developed for process engineers to quickly estimate and calculate and access a list of useful tables on site. As a part of a larger project to develop what is similar to Carl Branans book process engineers Pocket Handbook these tools would allow a process engineer to quickly calculate/estimate equipment efficiencies or sizing using standardized calculations The aim here is to take little-known data from the field (flow, pressures, temperatures, compositions..etc.) and use it as input for a rough estimation without having to return to the office to use commercial software or calculations Excel sheets to validate or to calculate. Additionally, these tools may also serve as a gathered data validation tool.'

colored_header.colored_header("Process Engineer's PocketApp",description,"red-70")

apps_list =[
"Heat exchanger sizing and rating",
"Physical properties",
"Line Sizing tool",
"Compressor Evaluation",
"Heater Efficiency",
"Cent. Pump performance test",
"Positive pump capacity calculations"

]

selected_app = st.selectbox('Select your calculation',apps_list, key = 'menu')

want_to_contribute = st.button("Let's Go!")
if want_to_contribute:
    switch_page(selected_app)