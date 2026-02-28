import streamlit as st
from streamlit_option_menu import option_menu
import base64
import re

def format_condition_label(label):
    """Format raw condition strings into user-friendly UI labels"""
    if not isinstance(label, str):
        label = str(label)
    
    # Replace underscores with spaces
    label = label.replace('_', ' ')
    
    # Strip leading weird artifacts like dots, bullet points, and whitespace
    label = re.sub(r'^[\.\s·\-]+', '', label)
    
    return label.strip()
def inject_custom_css():
    st.markdown("""
        <style>
        /* Hide the default Streamlit sidebar menu */
        [data-testid="stSidebarNav"] {
            display: none;
        }
        
        /* Main area background */
        .stApp {
            background-color: #f4f7f6;
        }
        
        /* Headers with gradient */
        h1, h2, h3 {
            color: #005b96;
            font-family: 'Inter', sans-serif;
            font-weight: 600;
        }
        h1 {
            background: linear-gradient(90deg, #005b96 0%, #0392cf 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            padding-bottom: 10px;
        }

        /* Styled Metric Cards */
        [data-testid="stMetricValue"] {
            color: #005b96;
            font-weight: 700;
            font-size: 2rem;
        }
        [data-testid="metric-container"] {
            background-color: white;
            border: 1px solid #e1e8ed;
            padding: 15px 20px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.04);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        [data-testid="metric-container"]:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 12px rgba(0, 91, 150, 0.1);
        }

        /* Buttons */
        .stButton>button {
            border-radius: 6px;
            font-weight: 500;
            transition: all 0.2s;
        }
        .stButton>button[kind="primary"] {
            background-color: #005b96;
            color: white;
        }
        .stButton>button[kind="primary"]:hover {
            background-color: #0392cf;
            border-color: #0392cf;
        }
        
        /* Dataframes */
        .stDataFrame {
            border: 1px solid #e1e8ed;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.02);
        }
        
        /* Text inputs & selects */
        .stTextInput>div>div>input {
            border-radius: 6px;
            border: 1px solid #cbd5e1;
        }
        .stTextInput>div>div>input:focus {
            border-color: #005b96;
            box-shadow: 0 0 0 1px #005b96;
        }
        
        /* Tabs */
        .stTabs [data-baseweb="tab-list"] {
            gap: 20px;
        }
        .stTabs [data-baseweb="tab"] {
            height: 50px;
            white-space: pre-wrap;
            border-radius: 6px 6px 0 0;
            padding: 10px 16px;
            color: #475569;
        }
        .stTabs [aria-selected="true"] {
            color: #005b96;
            background-color: #f0f7ff;
            border-bottom-color: #005b96;
        }
        
        /* Error/info banners */
        .stAlert {
            border-radius: 8px;
            border: none;
            box-shadow: 0 2px 5px rgba(0,0,0,0.03);
        }
        </style>
    """, unsafe_allow_html=True)

def render_sidebar_nav(current_page="Home"):
    inject_custom_css()
    
    with st.sidebar:
        st.markdown(
            """
            <div style='text-align: center; padding-bottom: 20px;'>
                <h2 style='color: #005b96; margin-bottom: 0;'><b>Midbase</b></h2>
                <div style='color: #64748b; font-size: 0.9em;'>Male Reproductive Genomics</div>
            </div>
            """, unsafe_allow_html=True
        )
        
        pages = ["Home", "Gene Search", "DE Explorer", "Upload & Compare", "Fetal Development", "Stem Cells", "Knowledge Graph", "About Us"]
        icons = ["house", "search", "bar-chart-line", "cloud-arrow-up", "activity", "diagram-3", "share", "info-circle"]
        
        try:
            default_ix = pages.index(current_page)
        except ValueError:
            default_ix = 0
            
        selected = option_menu(
            menu_title=None,
            options=pages,
            icons=icons,
            default_index=default_ix,
            styles={
                "container": {"padding": "0!important", "background-color": "transparent"},
                "icon": {"color": "#64748b", "font-size": "16px"}, 
                "nav-link": {"font-size": "14px", "text-align": "left", "margin":"2px", "color": "#1e2b3c", "border-radius": "6px", "padding": "10px"},
                "nav-link-selected": {"background-color": "#e0f0ff", "color": "#005b96", "font-weight": "600"},
            }
        )
        
        # Navigation logic using st.switch_page
        if selected != current_page:
            if selected == "Home":
                st.switch_page("Home.py")
            elif selected == "Gene Search":
                st.switch_page("pages/1_Gene_Search.py")
            elif selected == "DE Explorer":
                st.switch_page("pages/2_DE_Explorer.py")
            elif selected == "Upload & Compare":
                st.switch_page("pages/3_Upload_Compare.py")
            elif selected == "Fetal Development":
                st.switch_page("pages/4_Fetal_Development.py")
            elif selected == "Stem Cells":
                st.switch_page("pages/5_Stem_Cells.py")
            elif selected == "Knowledge Graph":
                st.switch_page("pages/6_Knowledge_Graph.py")
            elif selected == "About Us":
                st.switch_page("pages/7_About_Us.py")
