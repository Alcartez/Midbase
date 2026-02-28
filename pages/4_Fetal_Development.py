"""
Fetal Testis Development - GSE116278
Time-series gene expression across fetal development weeks 8-23
"""
import streamlit as st
import sqlite3
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

st.set_page_config(page_title="Fetal Development", page_icon="♂️", layout="wide")

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from style_utils import render_sidebar_nav
render_sidebar_nav("Fetal Development")

DB_PATH = "midbase.db"
PROCESSED_CSV = "data/processed/GSE116278.csv"

st.title("Fetal Testis Development")
st.markdown("**GSE116278** — Human fetal testis gene expression across developmental weeks **8–23**")

@st.cache_resource
def get_conn():
    return sqlite3.connect(DB_PATH, check_same_thread=False)

@st.cache_data
def load_fetal_data():
    df = pd.read_csv(PROCESSED_CSV, index_col=0, encoding='utf-8', encoding_errors='replace')
    return df

# Known key developmental genes
DEVELOPMENTAL_GENES = {
    "Sertoli Cell Markers": ["SOX9", "AMH", "FSHR", "GATA4", "WT1", "DHH"],
    "Leydig Cell Markers": ["STAR", "CYP11A1", "CYP17A1", "HSD3B2", "LHCGR"],
    "Germ Cell Markers": ["STRA8", "DAZL", "DDX4", "SYCP3", "PRDM9"],
    "Transcription Factors": ["NR5A1", "SOX17", "DMRT1", "PTGDS", "ARX"],
    "Signaling": ["WNT4", "FGF9", "GDNF", "INHBA"],
}

try:
    df = load_fetal_data()
    if 'Condition_Label' not in df.columns:
        st.error("Condition_Label column missing from data.")
        st.stop()
    
    # Extract week info from sample names / condition labels
    condition_col = df['Condition_Label']
    gene_df = df.drop(columns=['Condition_Label'])
    
    # Get unique time points (samples are rows)
    samples = df.index.tolist()
    conditions = condition_col.tolist()
    
    # Build a clean metadata table
    meta = pd.DataFrame({'sample': samples, 'condition': conditions})
    
    # Parse week number from sample name (e.g. "Fetal_Week_8_Day_0")
    def parse_week(sample_name):
        parts = str(sample_name).split('_')
        for i, p in enumerate(parts):
            if p == 'Week' and i+1 < len(parts):
                try:
                    return int(parts[i+1])
                except:
                    pass
        return None

    meta['week'] = meta['sample'].apply(parse_week)
    meta['original_idx'] = range(len(meta))
    meta = meta.dropna(subset=['week']).sort_values('week')
    week_order = sorted(meta['week'].unique())

    # Calculate gene variance across development weeks for dynamic selection
    # Filter for genes with mean > 0.5 to avoid zero-inflation noise
    overall_mean = gene_df.iloc[meta['original_idx']].mean()
    valid_genes = overall_mean[overall_mean > 0.5].index
    
    mean_by_week = gene_df[valid_genes].iloc[meta['original_idx']].copy()
    mean_by_week['week'] = meta['week'].values
    gene_variance = mean_by_week.groupby('week')[valid_genes].mean().std().sort_values(ascending=False)

    # ---- Sidebar ----
    with st.sidebar:
        st.header("Gene Selection")
        
        selection_mode = st.radio("Selection Method", ["Top Dynamic Genes", "Predefined Markers", "Custom"])
        
        if selection_mode == "Top Dynamic Genes":
            st.markdown("**Automatically select genes with the highest expression changes across development (Variance/SD):**")
            n_genes = st.select_slider("Select Top N dynamic genes", options=[5, 10, 25, 100], value=10)
            query_genes = gene_variance.head(n_genes).index.tolist()
            
        elif selection_mode == "Predefined Markers":
            st.markdown("**Quick select by cell type:**")
            selected_category = st.selectbox("Gene Category", list(DEVELOPMENTAL_GENES.keys()))
            query_genes = DEVELOPMENTAL_GENES[selected_category]
            
        else:
            st.markdown("**Enter custom gene(s):**")
            custom_input = st.text_area("Gene symbols (comma-separated)", "SOX9, AMH, FGF9")
            if custom_input.strip():
                query_genes = [g.strip().upper() for g in custom_input.split(',') if g.strip()]
            else:
                query_genes = ["SOX9"]
        
        st.markdown("---")
        plot_type = st.radio("Plot type", ["Line (mean per week)", "Box per week", "Heatmap"])
        normalize = st.checkbox("Normalize per gene (0–1)", value=False)

    # Verify genes exist in data
    available = [g for g in query_genes if g in gene_df.columns]
    missing = [g for g in query_genes if g not in gene_df.columns]

    if missing:
        st.warning(f"Not found in dataset: {', '.join(missing)}")
    if not available:
        st.error("None of the selected genes are in the dataset.")
        st.stop()

    # ---- Overview metrics ----
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Study", "GSE116278")
    col2.metric("Samples", len(meta))
    col3.metric("Time Points", len(week_order))
    col4.metric("Genes Available", f"{len(available)}/{len(query_genes)}")

    st.markdown("---")

    # ---- Build expression table per week ----
    # Join expression values with week metadata
    expr = gene_df[available].iloc[meta['original_idx']].copy()
    expr['week'] = meta['week'].values

    # ---- Plots ----
    if plot_type == "Heatmap":
        st.subheader("Expression Heatmap Across Development")
        
        mean_expr = expr.groupby('week')[available].mean()
        
        if normalize:
            # Min-max scale each gene
            mean_expr = (mean_expr - mean_expr.min()) / (mean_expr.max() - mean_expr.min() + 1e-9)
        
        fig = px.imshow(
            mean_expr.T,
            labels=dict(x="Week", y="Gene", color="Expression"),
            color_continuous_scale="RdBu_r",
            aspect="auto",
            title="Mean Expression per Developmental Week",
            x=[f"Wk {w}" for w in mean_expr.index],
        )
        fig.update_layout(height=max(300, len(available) * 40))
        st.plotly_chart(fig, width='stretch')

    elif plot_type == "Line (mean per week)":
        st.subheader("Expression Trajectories Across Development")
        
        mean_expr = expr.groupby('week')[available].mean().reset_index()
        
        melted = mean_expr.melt(id_vars='week', var_name='Gene', value_name='Expression')
        
        if normalize:
            melted['Expression'] = melted.groupby('Gene')['Expression'].transform(
                lambda x: (x - x.min()) / (x.max() - x.min() + 1e-9)
            )
        
        fig = px.line(
            melted, x='week', y='Expression', color='Gene',
            markers=True,
            title="Gene Expression Across Fetal Development",
            labels={'week': 'Gestational Week', 'Expression': 'Mean Expression' + (' (normalized)' if normalize else '')},
        )
        fig.update_layout(height=500, xaxis=dict(tickmode='linear', dtick=1))
        st.plotly_chart(fig, width='stretch')

    else:  # Box per week
        st.subheader("Expression Distribution by Week")
        
        gene_choice = st.selectbox("Select gene for boxplot", available)
        plot_data = expr[['week', gene_choice]].rename(columns={gene_choice: 'Expression'})
        
        fig = px.box(
            plot_data, x='week', y='Expression',
            title=f"{gene_choice} — Expression by Gestational Week",
            labels={'week': 'Gestational Week'},
            color='week', color_continuous_scale='Blues',
        )
        fig.update_layout(height=450, showlegend=False)
        st.plotly_chart(fig, width='stretch')

    st.markdown("---")

    # ---- Top dynamic genes (most variable across weeks) ----
    st.subheader("Most Dynamically Expressed Genes Across Development")
    top_dynamic = gene_variance.head(20).reset_index()
    top_dynamic.columns = ['Gene', 'Expression SD across weeks']
    top_dynamic['Expression SD across weeks'] = top_dynamic['Expression SD across weeks'].round(3)

    col_a, col_b = st.columns([1, 2])
    with col_a:
        st.dataframe(top_dynamic, hide_index=True, width='stretch')
    with col_b:
        fig2 = px.bar(
            top_dynamic, x='Gene', y='Expression SD across weeks',
            title="Top 20 Most Variable Genes Across Development",
            color='Expression SD across weeks', color_continuous_scale='Blues'
        )
        fig2.update_layout(height=400)
        st.plotly_chart(fig2, width='stretch')

    # ---- Download ----
    st.markdown("---")
    mean_expr_download = mean_by_week.groupby('week')[valid_genes].mean()
    csv = mean_expr_download.to_csv()
    st.download_button("Download Mean Expression by Week (CSV)", csv,
                       file_name="GSE116278_expression_by_week.csv", mime="text/csv")

except FileNotFoundError:
    st.error(f"Data file not found: {PROCESSED_CSV}")
except Exception as e:
    st.error(f"Error loading data: {e}")
    import traceback; st.code(traceback.format_exc())
