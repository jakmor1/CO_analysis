import streamlit as st

st.title("Application for analyzing the&nbsp;distribution&nbsp;of crossing-over in&nbsp;Arabidopsis thaliana")

st.header("General information")
st.write("Welcome to the webpage dedicated to the analysis of crossing-over distribution in Arabidopsis thaliana. This project was undertaken as part of a bachelor's thesis, conducted under the guidance of Professor Piotr Ziółkowski, and with the invaluable support of Wojciech Dzięgielewski. The data used for the analysis was obtained using the Trained Individual GenomE Reconstruction (Rowan et al. 2015).")

st.header("Contact")
st.write("Email correspondence can be sent to:  \njakub.morawski@student.put.poznan.pl, pzio@amu.edu.pl and wojdzi@amu.edu.pl ")


# Stopka
st.write("---")
st.write("©2023 Institute of Molecular Biology and Biotechnology, Laboratory of Genome Biology, Adam&nbsp;Mickiewicz&nbsp;University, Poznan, Poland ")


hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 
