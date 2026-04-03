"""
Differential Expression Explorer - Volcano plots and heatmaps
Updated to handle multiple comparisons per study and fix gene count discrepancies.
"""
import streamlit as st
import sqlite3
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from pathlib import Path

st.set_page_config(page_title="DE Explorer", page_icon="♂️", layout="wide")

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from style_utils import render_sidebar_nav
render_sidebar_nav("DE Explorer")

st.title("Differential Expression Explorer")

# Database connections
@st.cache_resource
def get_de_db_connection():
    db_path = Path("midbase_de.db")
    return sqlite3.connect(str(db_path), check_same_thread=False)

@st.cache_resource
def get_core_db_connection():
    db_path = Path("midbase_core.db")
    return sqlite3.connect(str(db_path), check_same_thread=False)

@st.cache_data(ttl=3600)
def get_studies():
    conn = get_core_db_connection()
    # Exclude special dataset pages - Fetal Timeseries and Stem Cells have their own pages
    EXCLUDED = ('GSE116278', 'GSE92276')
    df = pd.read_sql_query("SELECT study_id, title FROM studies ORDER BY study_id", conn)
    return df[~df['study_id'].isin(EXCLUDED)]

@st.cache_data(ttl=3600)
def get_comparisons(study_id):
    """Fetch unique comparisons for a specific study to prevent duplicate gene counts."""
    conn = get_de_db_connection()
    query = f"SELECT DISTINCT comparison FROM differential_expression WHERE study_id = '{study_id}'"
    df = pd.read_sql_query(query, conn)
    return df['comparison'].tolist()

@st.cache_data(ttl=3600)
def get_de_results(study_id, comparison, p_threshold=1.0):
    """Fetch DE results filtered by both study and the specific comparison group."""
    conn = get_de_db_connection()
    sql = f"""
    SELECT gene_symbol, logFC, ave_expr, p_value, comparison
    FROM differential_expression
    WHERE study_id = '{study_id}' 
    AND comparison = ?
    AND p_value <= {p_threshold}
    ORDER BY p_value
    """
    return pd.read_sql_query(sql, conn, params=(comparison,))

# Sidebar filters
with st.sidebar:
    st.header("Filters")
    
    studies = get_studies()
    selected_study = st.selectbox(
        "Select Study",
        studies['study_id'].tolist(),
        format_func=lambda x: f"{x} - {studies[studies['study_id']==x]['title'].iloc[0][:40]}"
    )

    # Comparison Selector - Fixes the 67k vs 134k gene count issue
    available_comparisons = get_comparisons(selected_study)
    if available_comparisons:
        selected_comparison = st.selectbox("Select Comparison Group", available_comparisons)
    else:
        selected_comparison = None
        st.warning("No specific comparisons found for this study.")
    
    st.markdown("---")
    
    p_threshold = st.slider(
        "P-value threshold",
        min_value=0.0,
        max_value=1.0,
        value=1.0,
        step=0.01
    )
    
    logfc_threshold = st.slider(
        "|LogFC| threshold",
        min_value=0.0,
        max_value=5.0,
        value=0.0,
        step=0.1
    )
    
    st.markdown("---")
    show_labels = st.checkbox("Show gene labels (top 10)", value=False)

# Load data based on selection
if selected_comparison:
    de_data = get_de_results(selected_study, selected_comparison, p_threshold)
else:
    de_data = pd.DataFrame()

if len(de_data) == 0:
    st.warning(f"No differential expression results found for {selected_study} ({selected_comparison})")
else:
    # Apply logFC filter
    de_data_filtered = de_data[de_data['logFC'].abs() >= logfc_threshold].copy()
    
    if len(de_data_filtered) == 0:
        st.warning(f"No genes pass the logFC threshold of {logfc_threshold}")
    else:
        # Stats
        st.subheader(f"{selected_study} - {selected_comparison}")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Genes", len(de_data))
        with col2:
            sig_count = (de_data['p_value'] < 0.05).sum()
            st.metric("Significant (p<0.05)", sig_count)
        with col3:
            up_count = ((de_data['logFC'] > 1) & (de_data['p_value'] < 0.05)).sum()
            st.metric("Upregulated", up_count)
        with col4:
            down_count = ((de_data['logFC'] < -1) & (de_data['p_value'] < 0.05)).sum()
            st.metric("Downregulated", down_count)
        
        st.markdown("---")
        
        # Volcano plot
        st.subheader("Volcano Plot")
        
        # Color by significance
        de_data_filtered['Significant'] = 'No'
        de_data_filtered.loc[(de_data_filtered['p_value'] < 0.05) & (de_data_filtered['logFC'].abs() > 1), 'Significant'] = 'Yes'
        
        fig = px.scatter(
            de_data_filtered,
            x='logFC',
            y='p_value',
            color='Significant',
            color_discrete_map={'Yes': 'red', 'No': 'gray'},
            hover_name='gene_symbol',
            hover_data={'p_value': ':.2e', 'logFC': ':.2f', 'ave_expr': ':.2f'},
            opacity=0.6,
            title=f"Volcano Plot: {selected_study} ({selected_comparison})",
            labels={'logFC': 'Log2 Fold Change', 'p_value': 'P-value'}
        )
        
        # Add threshold lines
        fig.add_hline(y=0.05, line_dash="dash", line_color="black", annotation_text="p=0.05")
        fig.add_vline(x=1, line_dash="dash", line_color="black")
        fig.add_vline(x=-1, line_dash="dash", line_color="black")
        
        # Add labels for top genes
        if show_labels:
            top_genes = de_data_filtered.nsmallest(10, 'p_value')
            for _, row in top_genes.iterrows():
                fig.add_annotation(
                    x=row['logFC'],
                    y=row['p_value'],
                    text=row['gene_symbol'],
                    showarrow=True,
                    arrowhead=2,
                    arrowsize=1,
                    arrowwidth=1,
                    arrowcolor="#636363",
                    ax=20,
                    ay=-30
                )
        
        fig.update_layout(height=600)
        st.plotly_chart(fig, use_container_width=True)
        
        st.markdown("---")
        
        # Top genes table
        st.subheader("Top Differentially Expressed Genes")
        
        top_n = st.slider("Number of genes to display", 10, 100, 20)
        
        top_genes_df = de_data_filtered.nsmallest(top_n, 'p_value')[['gene_symbol', 'logFC', 'ave_expr', 'p_value', 'comparison']]
        
        # Format for display
        top_genes_display = top_genes_df.copy()
        top_genes_display['p_value'] = top_genes_display['p_value'].apply(lambda x: f"{x:.2e}")
        top_genes_display['logFC'] = top_genes_display['logFC'].apply(lambda x: f"{x:.3f}")
        top_genes_display['ave_expr'] = top_genes_display['ave_expr'].apply(lambda x: f"{x:.2f}")
        
        st.dataframe(top_genes_display, use_container_width=True)
        
        # Download
        csv = top_genes_df.to_csv(index=False)
        st.download_button(
            label=f"Download Top {top_n} Genes (CSV)",
            data=csv,
            file_name=f"{selected_study}_{selected_comparison}_top_{top_n}_genes.csv",
            mime="text/csv"
        )
