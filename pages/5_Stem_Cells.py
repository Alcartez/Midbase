"""
Spermatogonial Stem Cells - GSE92276
175 sorted SSC samples — marker gene expression and cell identity
"""
import streamlit as st
import sqlite3
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

st.set_page_config(page_title="Stem Cells", page_icon="♂️", layout="wide")

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from style_utils import render_sidebar_nav
render_sidebar_nav("Stem Cells")

DB_PATH  = "midbase.db"
PROCESSED_CSV = "data/processed/GSE92276.csv"

st.title("Spermatogonial Stem Cells")
st.markdown("**GSE92276** — 175 sorted human spermatogonial stem cell (SSC) samples")

# Known SSC biology gene sets
SSC_GENE_SETS = {
    "SSC Self-Renewal": ["GDNF", "GFRA1", "RET", "NANOS3", "ID4", "BCAS3"],
    "SSC Markers": ["DAZL", "DDX4", "STRA8", "UTF1", "UCHL1", "THY1"],
    "Differentiation": ["KIT", "STRA8", "SYCP1", "SYCP3", "PRDM9", "MEIOB"],
    "Transcription Factors": ["DMRT1", "NANOS2", "SOHLH1", "SOHLH2", "ZBTB16", "ID4"],
    "Signaling": ["GDNF", "FGF2", "BMP4", "SMAD1", "SMAD5"],
    "Epigenetic": ["DNMT3A", "DNMT3L", "PIWIL4", "UHRF1", "SETDB1"],
}

@st.cache_data
def load_ssc_data():
    df = pd.read_csv(PROCESSED_CSV, index_col=0, encoding='utf-8', encoding_errors='replace')
    return df

try:
    df = load_ssc_data()
    if 'Condition_Label' in df.columns:
        gene_df = df.drop(columns=['Condition_Label'])
    else:
        gene_df = df

    # Pre-compute metrics to avoid zero-inflation noise (filter mean > 0.5)
    overall_mean_all = gene_df.mean()
    valid_genes = overall_mean_all[overall_mean_all > 0.5].index
    gene_cv_all = (gene_df[valid_genes].std() / gene_df[valid_genes].mean()).sort_values(ascending=False)
    overall_mean_sorted = overall_mean_all.sort_values(ascending=False)

    # ---- Sidebar ----
    with st.sidebar:
        st.header("Gene Selection")

        selection_mode = st.radio("Selection Method", [
            "Top Expressed Genes", 
            "Predefined Markers", 
            "Custom"
        ])
        
        if selection_mode == "Top Expressed Genes":
            st.markdown("**Select genes with the highest overall mean expression:**")
            n_genes = st.select_slider("Select Top N genes", options=[5, 10, 25, 100], value=5)
            overall_mean = overall_mean_sorted.head(n_genes)
            query_genes = overall_mean.index.tolist()
            selected_set = f"Top {n_genes} Expressed Genes"
            

            
        elif selection_mode == "Predefined Markers":
            st.markdown("**Quick select by known biology:**")
            selected_set = st.selectbox("Gene Set", list(SSC_GENE_SETS.keys()))
            query_genes = SSC_GENE_SETS[selected_set]
            
        else:
            st.markdown("**Enter custom gene(s):**")
            custom = st.text_area("Gene symbols (comma-separated)", "GDNF, GFRA1, ZBTB16")
            selected_set = "Custom Gene List"
            query_genes = [g.strip().upper() for g in custom.split(',') if g.strip()] if custom.strip() else ["GDNF"]

        st.markdown("---")
        st.markdown("**Visualization**")
        n_top = st.slider("Top expressed genes to show", 10, 100, 30)

    available = [g for g in query_genes if g in gene_df.columns]
    missing   = [g for g in query_genes if g not in gene_df.columns]

    if missing:
        st.warning(f"Not found: {', '.join(missing)}")

    # ---- Overview metrics ----
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Study", "GSE92276")
    col2.metric("Samples", len(df))
    col3.metric("Total Genes", len(gene_df.columns))
    col4.metric("Genes Matched", f"{len(available)}/{len(query_genes)}")

    st.markdown("---")

    # ---- Mean expression of selected gene set ----
    if available:
        st.subheader(f"{selected_set} — Expression Across 175 SSC Samples")

        expr_subset = gene_df[available]
        mean_vals = expr_subset.mean().sort_values(ascending=False).reset_index()
        mean_vals.columns = ['Gene', 'Mean Expression']

        col_left, col_right = st.columns([1, 2])
        with col_left:
            display = mean_vals.copy()
            display['Mean Expression'] = display['Mean Expression'].round(2)
            st.dataframe(display, hide_index=True, width='stretch')

        with col_right:
            fig = px.bar(
                mean_vals, x='Gene', y='Mean Expression',
                color='Mean Expression', color_continuous_scale='Blues',
                title=f"Mean Expression: {selected_set}",
            )
            fig.update_layout(height=400)
            st.plotly_chart(fig, width='stretch')

        st.markdown("---")

        # ---- Distribution violin/box per gene ----
        st.subheader("Expression Distribution Across Samples")

        melted = expr_subset.melt(var_name='Gene', value_name='Expression')
        fig2 = px.violin(
            melted, x='Gene', y='Expression',
            box=True, points=False,
            color='Gene',
            title=f"Expression Distribution: {selected_set}",
        )
        fig2.update_layout(height=450, showlegend=False)
        st.plotly_chart(fig2, width='stretch')

    st.markdown("---")

    # ---- Top expressed genes overall ----
    st.subheader(f"Top {n_top} Highest Expressed Genes Across All SSC Samples")

    overall_mean = overall_mean_sorted.head(n_top).reset_index()
    overall_mean.columns = ['Gene', 'Mean Expression']
    overall_mean['Mean Expression'] = overall_mean['Mean Expression'].round(2)

    col_a, col_b = st.columns([1, 2])
    with col_a:
        st.dataframe(overall_mean, hide_index=True, width='stretch')
    with col_b:
        fig3 = px.bar(
            overall_mean, x='Gene', y='Mean Expression',
            color='Mean Expression', color_continuous_scale='Blues',
            title=f"Top {n_top} Expressed Genes in SSCs",
        )
        fig3.update_layout(height=450)
        st.plotly_chart(fig3, width='stretch')



    st.markdown("---")

    # ---- Download ----
    summary = gene_df.agg(['mean', 'std', 'median']).T
    summary.columns = ['Mean', 'Std', 'Median']
    summary = summary.round(3)
    csv = summary.to_csv()
    st.download_button("Download SSC Expression Summary (CSV)", csv,
                       file_name="GSE92276_SSC_expression_summary.csv", mime="text/csv")

except FileNotFoundError:
    st.error(f"Data file not found: {PROCESSED_CSV}")
except Exception as e:
    st.error(f"Error loading data: {e}")
    import traceback; st.code(traceback.format_exc())
