"""
XY Counsel - Testicular Gene Expression Database Platform
Home Page: Dashboard and overview statistics
"""
import streamlit as st
import sqlite3
import pandas as pd
from pathlib import Path

# Page config
st.set_page_config(
    page_title="Midbase - Gene Expression Database",
    page_icon="♂️",
    layout="wide",
    initial_sidebar_state="expanded"
)

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))
from style_utils import render_sidebar_nav

render_sidebar_nav("Home")

# Custom CSS for Home
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        background: linear-gradient(90deg, #005b96 0%, #0392cf 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 1rem;
    }
    .stat-card {
        background: linear-gradient(135deg, #005b96 0%, #0392cf 100%);
        padding: 1.5rem;
        border-radius: 10px;
        color: white;
        text-align: center;
        box-shadow: 0 4px 6px rgba(0, 91, 150, 0.2);
    }
    .stat-number {
        font-size: 2.5rem;
        font-weight: bold;
    }
    .stat-label {
        font-size: 1rem;
        opacity: 0.9;
    }
    .study-card {
        border: 1px solid #e1e8ed;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
        background: white;
        transition: transform 0.2s, box-shadow 0.2s;
    }
    .study-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 8px rgba(0, 91, 150, 0.1);
    }
    .study-title {
        font-weight: bold;
        color: #005b96;
    }
</style>
""", unsafe_allow_html=True)

# Database connection
@st.cache_resource
def get_database_connection():
    db_path = Path("midbase_core.db")
    if not db_path.exists():
        st.error(f"Database not found: {db_path} (Are you running from the right directory?)")
        return None
    return sqlite3.connect(str(db_path), check_same_thread=False)

@st.cache_data(ttl=3600)
def get_database_stats():
    """Get overview statistics from database"""
    conn = get_database_connection()
    if conn is None:
        return None
    
    stats = {}
    
    # Number of studies
    stats['studies'] = pd.read_sql_query("SELECT COUNT(*) as count FROM studies", conn).iloc[0]['count']
    
    # Number of samples
    stats['samples'] = pd.read_sql_query("SELECT COUNT(*) as count FROM samples", conn).iloc[0]['count']
    
    # Connect to the split DE database for the remaining stats
    try:
        conn_de = sqlite3.connect("midbase_de.db")
        # Number of DE results
        stats['de_results'] = pd.read_sql_query("SELECT COUNT(*) as count FROM differential_expression", conn_de).iloc[0]['count']
        
        # Unique genes
        stats['genes'] = pd.read_sql_query("SELECT COUNT(DISTINCT gene_symbol) as count FROM differential_expression", conn_de).iloc[0]['count']
        conn_de.close()
    except Exception as e:
        stats['de_results'] = 0
        stats['genes'] = 0
        st.warning(f"Could not load DE stats: {e}")
        
    return stats

@st.cache_data(ttl=3600)
def get_all_studies():
    """Get all studies with metadata"""
    conn = get_database_connection()
    if conn is None:
        return None
    
    query = """
    SELECT 
        s.study_id,
        s.title,
        s.category,
        COUNT(DISTINCT sa.sample_id) as sample_count
    FROM studies s
    LEFT JOIN samples sa ON s.study_id = sa.study_id
    GROUP BY s.study_id, s.title, s.category
    ORDER BY s.category, s.study_id
    """
    return pd.read_sql_query(query, conn)

# Header
st.markdown('<p class="main-header">Midbase (Project XY_Counsel)</p>', unsafe_allow_html=True)
st.markdown("### Testicular Gene Expression Database & Analysis Platform")

st.markdown("---")

# Overview statistics
st.subheader("Database Overview")

stats = get_database_stats()
if stats:
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-number">{stats['studies']}</div>
            <div class="stat-label">Studies</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-number">{stats['samples']}</div>
            <div class="stat-label">Samples</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-number">{stats['genes']:,}</div>
            <div class="stat-label">Genes Analyzed</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-number">{stats['de_results']:,}</div>
            <div class="stat-label">DE Results</div>
        </div>
        """, unsafe_allow_html=True)

st.markdown("---")

# Study browser
st.subheader("📚 Available Datasets")

studies_df = get_all_studies()
if studies_df is not None and len(studies_df) > 0:
    # Group by category
    categories = studies_df['category'].unique()
    
    for category in sorted(categories):
        st.markdown(f"**{category}**")
        category_studies = studies_df[studies_df['category'] == category]
        
        for _, study in category_studies.iterrows():
            with st.container():
                st.markdown(f"""
                <div class="study-card">
                    <div class="study-title">{study['study_id']}</div>
                    <div>{study['title']}</div>
                    <div style="color: #666; font-size: 0.9rem; margin-top: 0.5rem;">
                        {study['sample_count']} samples
                    </div>
                </div>
                """, unsafe_allow_html=True)
        
        st.markdown("")

else:
    st.warning("No studies found in database.")

# Quick search
st.markdown("---")
st.subheader("Quick Gene Search")

col_search, col_button = st.columns([4, 1])
with col_search:
    gene_query = st.text_input("Enter gene symbol (e.g., TP53, BRCA1)", placeholder="Gene symbol...")
with col_button:
    st.markdown("<br>", unsafe_allow_html=True)
    search_clicked = st.button("Search", use_container_width=True)

if search_clicked and gene_query:
    conn = get_database_connection()
    if conn:
        query = f"""
        SELECT study_id, gene_symbol, logFC, p_value, comparison
        FROM differential_expression
        WHERE gene_symbol LIKE '%{gene_query}%'
        ORDER BY p_value
        LIMIT 10
        """
        results = pd.read_sql_query(query, conn)
        
        if len(results) > 0:
            st.success(f"Found {len(results)} results for '{gene_query}'")
            st.dataframe(results, use_container_width=True)
        else:
            st.info(f"No results found for '{gene_query}'")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; font-size: 0.9rem;">
    <p>XY Counsel Platform | 9 Datasets | 301 Samples | 67,186 Genes</p>
    <p>Use the sidebar to navigate between pages</p>
</div>
""", unsafe_allow_html=True)
