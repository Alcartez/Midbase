import streamlit as st
import sys
from pathlib import Path

st.set_page_config(page_title="About Us", page_icon="ℹ️", layout="wide")

sys.path.append(str(Path(__file__).parent.parent))
from style_utils import render_sidebar_nav

render_sidebar_nav("About Us")

st.title("About Us")
st.markdown("---")

col1, col2 = st.columns([2, 1])

with col1:
    st.header("SBDA Research Lab")
    st.subheader("Amity University, Noida")
    
    st.markdown("""
    Welcome to Midbase, a comprehensive platform for male reproductive genomics developed at the **Systems Biology and Data Analytics (SBDA) Research Lab** at **Amity University, Noida**.
    
    Our lab is dedicated to advancing the understanding of complex biological systems through computational analysis, bioinformatics, and machine learning. Midbase represents our commitment to democratizing access to high-quality transcriptomic data and analytical tools for researchers worldwide.
    
    ### Our Mission
    - **Data Integration:** Bridging the gap between disparate transcriptomic datasets from global repositories.
    - **Advanced Diagnostics:** Developing computational pipelines (like our Diagnostic Knowledge Graph and PCA Comparative Analysis) to identify molecular markers for male infertility and reproductive conditions.
    - **Accessibility:** Providing a unified, user-friendly computational interface so that biologists and clinicians can perform complex analyses without needing to write code.
    
    ### Research Focus Areas
    - Male Infertility & Reproductive Health
    - Bioinformatics Pipeline Development
    - AI / ML in Genomics
    - Systems Biology and Regulatory Networks
    """)
    
with col2:
    st.info("""
    **Contact Information**
    
    **Institution:** 
    Amity University Uttar Pradesh
    Sector 125, Noida 201313
    India
    
    **Lab:**
    Systems Biology and Data Analytics (SBDA) Research Lab
    """)
    
st.markdown("---")
st.markdown("<p style='text-align: center; color: gray;'>Designed and developed with 💻 & 🧬 by the SBDA Lab</p>", unsafe_allow_html=True)
