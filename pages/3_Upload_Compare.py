"""
Upload & Compare - User Data Analysis
Upload your own gene expression data and compare against midbase database
"""
import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import os
import plotly.express as px
import plotly.graph_objects as go
import io

# Manual implementations (no sklearn/scipy due to Windows DLL blocking)
def manual_pca(X, n_components=2):
    """Manual PCA implementation using SVD (memory efficient)"""
    # Center the data
    X_centered = X - np.mean(X, axis=0)
    
    # Use SVD instead of covariance matrix (much more memory efficient!)
    # For X with shape (n_samples, n_features)
    # SVD gives us: X = U @ S @ Vt
    # Principal components are in Vt
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    
    # Project data onto principal components
    # X_pca = X_centered @ Vt.T[:, :n_components]
    # But we can use U @ diag(S) which is already computed
    X_pca = U[:, :n_components] * S[:n_components]
    
    # Calculate explained variance ratio
    explained_variance = (S ** 2) / (X_centered.shape[0] - 1)
    total_var = np.sum(explained_variance)
    explained_variance_ratio = explained_variance[:n_components] / total_var
    
    return X_pca, explained_variance_ratio

def standardize(X):
    """Standardize data (z-score normalization)"""
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)
    std[std == 0] = 1  # Avoid division by zero
    return (X - mean) / std

def pearson_correlation(x, y):
    """Calculate Pearson correlation coefficient"""
    x_centered = x - np.mean(x)
    y_centered = y - np.mean(y)
    
    numerator = np.sum(x_centered * y_centered)
    denominator = np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2))
    
    if denominator == 0:
        return 0.0
    
    return numerator / denominator

# Page config
st.set_page_config(page_title="Upload & Compare", page_icon="♂️", layout="wide")

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from style_utils import render_sidebar_nav, format_condition_label
render_sidebar_nav("Upload & Compare")

# Custom CSS (overrides globals if necessary, but keep base gradient if preferred)
st.markdown("""
<style>
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
</style>
""", unsafe_allow_html=True)

st.title("Upload & Compare Your Data")
st.markdown("Upload your gene expression data to compare against our testicular gene expression database")

# Database path
DB_PATH = "midbase_core.db"
PROCESSED_DIR = "data/processed"

@st.cache_data
def load_midbase_expression():
    """Load all expression data from processed CSVs"""
    all_data = []
    metadata = []
    
    # Get list of all processed files
    csv_files = [f for f in os.listdir(PROCESSED_DIR) if f.endswith('.csv')]
    
    for csv_file in csv_files:
        study_id = csv_file.replace('.csv', '')
        csv_path = os.path.join(PROCESSED_DIR, csv_file)
        
        try:
            df = pd.read_csv(csv_path, index_col=0)
            
            # Check if has Condition_Label
            if 'Condition_Label' in df.columns:
                conditions = df['Condition_Label']
                df = df.drop(columns=['Condition_Label'])
                
                # Store metadata
                for sample_name in df.index:
                    raw_cond = conditions.loc[sample_name]
                    # Handle duplicated index bugs returning pandas Series objects
                    if isinstance(raw_cond, pd.Series):
                        raw_cond = str(raw_cond.iloc[0])
                    else:
                        raw_cond = str(raw_cond)
                        
                    metadata.append({
                        'sample_id': sample_name,
                        'study_id': study_id,
                        'condition': format_condition_label(raw_cond)
                    })
            else:
                # No condition labels
                for sample_name in df.index:
                    metadata.append({
                        'sample_id': sample_name,
                        'study_id': study_id,
                        'condition': 'Unknown'
                    })
            
            all_data.append(df)
        except Exception as e:
            st.warning(f"Skipping {study_id}: {e}")
            continue
    
    # Concatenate vertically (samples as rows, genes as columns)
    combined = pd.concat(all_data, axis=0)
    metadata_df = pd.DataFrame(metadata)
    
    return combined, metadata_df

def validate_user_data(df):
    """Validate uploaded data"""
    errors = []
    warnings = []
    
    # Check minimum genes
    if len(df.index) < 1000:
        errors.append(f"Too few genes ({len(df.index)}). Minimum 1,000 genes required for reliable comparison.")
    
    # Check for numeric data
    if not df.select_dtypes(include=[np.number]).shape[1] > 0:
        errors.append("No numeric columns found. Data must contain expression values.")
    
    # Check for gene symbols in index
    if df.index.dtype != 'object':
        warnings.append("Gene names should be in the first column (index).")
    
    # Check data range (detect if log-transformed)
    mean_val = df.mean().mean()
    if mean_val > 100:
        warnings.append("Data appears to be raw counts. Consider log2 transformation.")
    elif mean_val < 0:
        errors.append("Negative values detected. Expression data should be non-negative.")
    
    return {'valid': len(errors) == 0, 'errors': errors, 'warnings': warnings}

def preprocess_user_data(df, auto_log=False):
    """Preprocess user uploaded data"""
    # Drop non-numeric columns
    df_numeric = df.select_dtypes(include=[np.number])
    
    # Auto log-transform if requested and data looks like raw counts
    if auto_log and df_numeric.mean().mean() > 50:
        st.info("Applying log2(x+1) transformation...")
        df_numeric = np.log2(df_numeric + 1)
    
    # Expected upload format: genes as rows, samples as columns
    # We need: samples as rows, genes as columns (to match midbase)
    # So always transpose
    if len(df_numeric.index) > len(df_numeric.columns):
        # Likely genes are rows (many genes, few samples)
        st.info(f"Transposing data: {len(df_numeric.index)} genes × {len(df_numeric.columns)} samples → {len(df_numeric.columns)} samples × {len(df_numeric.index)} genes")
        df_numeric = df_numeric.T
    else:
        # Already in correct format (samples as rows)
        st.info(f"Data format OK: {len(df_numeric.index)} samples × {len(df_numeric.columns)} genes")
    
    return df_numeric

def calculate_similarity(user_data, midbase_data, metadata_df, top_n=10):
    """Calculate correlation-based similarity scores"""
    # Exclude conditions that distort disease correlations
    exclude_conditions = ['Fetal Timeseries', 'Spermatogonial Stem Cells']
    filtered_meta = metadata_df[~metadata_df['condition'].isin(exclude_conditions)]
    filtered_samples = filtered_meta['sample_id'].tolist()
    
    valid_midbase_samples = [s for s in midbase_data.index if s in filtered_samples]
    filtered_midbase = midbase_data.loc[valid_midbase_samples]
    
    # Find common genes
    common_genes = user_data.columns.intersection(filtered_midbase.columns)
    
    if len(common_genes) < 500:
        st.error(f"Too few common genes ({len(common_genes)}). Need at least 500 for comparison.")
        return None
    
    st.success(f"Found {len(common_genes):,} common genes for comparison")
    
    # Subset to common genes
    user_common = user_data[common_genes]
    midbase_common = filtered_midbase[common_genes]
    
    results = []
    
    with st.spinner("Calculating similarity scores..."):
        for user_sample in user_common.index:
            user_vec = user_common.loc[user_sample].values
            
            correlations = []
            for midbase_sample in midbase_common.index:
                midbase_vec = midbase_common.loc[midbase_sample].values
                
                # Calculate Pearson correlation (manual)
                corr = pearson_correlation(user_vec, midbase_vec)
                
                # Get metadata
                meta = metadata_df[metadata_df['sample_id'] == midbase_sample].iloc[0]
                
                correlations.append({
                    'user_sample': user_sample,
                    'midbase_sample': midbase_sample,
                    'study_id': meta['study_id'],
                    'condition': meta['condition'],
                    'correlation': corr
                })
            
            # Sort by correlation and keep top N
            correlations = sorted(correlations, key=lambda x: x['correlation'], reverse=True)
            results.extend(correlations[:top_n])
    
    return pd.DataFrame(results)

def plot_pca_overlay(midbase_data, user_data, metadata_df):
    """PCA with user data overlaid (manual implementation)"""
    # Exclude conditions that distort disease clustering
    exclude_conditions = ['Fetal Timeseries', 'Spermatogonial Stem Cells']
    filtered_meta = metadata_df[~metadata_df['condition'].isin(exclude_conditions)]
    filtered_samples = filtered_meta['sample_id'].tolist()
    
    # Keep only samples not in exclude list
    valid_midbase_samples = [s for s in midbase_data.index if s in filtered_samples]
    filtered_midbase = midbase_data.loc[valid_midbase_samples]
    
    # Find common genes
    common_genes = user_data.columns.intersection(filtered_midbase.columns)
    
    # Combine data (samples as rows, genes as columns)
    combined = pd.concat([
        filtered_midbase[common_genes],
        user_data[common_genes]
    ], axis=0)
    
    # Standardize
    X_scaled = standardize(combined.values)
    
    # Manual PCA
    X_pca, explained_var = manual_pca(X_scaled, n_components=2)
    
    # Create labels
    labels = []
    colors = []
    sizes = []
    symbols = []
    for idx in combined.index:
        if idx in user_data.index:
            labels.append(f'Your Data: {idx}')
            colors.append('User Data')
            sizes.append(15)
            symbols.append('User Upload')
        else:
            meta = metadata_df[metadata_df['sample_id'] == idx]
            if len(meta) > 0:
                condition = meta.iloc[0]['condition']
                labels.append(f"{meta.iloc[0]['study_id']}: {condition}")
                colors.append(condition)
            else:
                labels.append(idx)
                colors.append('Unknown')
            sizes.append(8)
            symbols.append('Reference Database')
    
    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1],
        'Sample': labels,
        'Type': colors,
        'Size': sizes,
        'Source': symbols
    })
    
    # Plot
    fig = px.scatter(
        plot_df,
        x='PC1',
        y='PC2',
        color='Type',
        color_discrete_sequence=px.colors.qualitative.Alphabet,
        symbol='Source',
        symbol_map={'User Upload': 'star', 'Reference Database': 'circle'},
        hover_name='Sample',
        size='Size',
        size_max=15,
        title=f'PCA: Your Data vs Midbase (PC1: {explained_var[0]:.1%}, PC2: {explained_var[1]:.1%})',
        labels={'PC1': f'PC1 ({explained_var[0]:.1%})',
                'PC2': f'PC2 ({explained_var[1]:.1%})'},
        height=600
    )
    
    fig.update_traces(marker=dict(line=dict(width=1, color='white')))
    
    return fig

def plot_correlation_heatmap(similarity_df, user_samples):
    """Plot correlation heatmap showing top matches"""
    # Pivot data for heatmap
    pivot_data = []
    
    for user_sample in user_samples:
        sample_data = similarity_df[similarity_df['user_sample'] == user_sample].head(20)
        for _, row in sample_data.iterrows():
            pivot_data.append({
                'User Sample': user_sample,
                'Midbase Sample': f"{row['study_id']}: {row['condition']}",
                'Correlation': row['correlation']
            })
    
    pivot_df = pd.DataFrame(pivot_data)
    
    if len(pivot_df) == 0:
        return None
    
    # Create heatmap
    fig = px.density_heatmap(
        pivot_df,
        x='Midbase Sample',
        y='User Sample',
        z='Correlation',
        color_continuous_scale='Blues',
        title='Top 20 Matches: Correlation Heatmap',
        height=400
    )
    
    fig.update_xaxes(tickangle=45)
    
    return fig

def calculate_kg_projection(user_data, kg_path="data/knowledge_graph.json"):
    """Project user samples onto the Diagnostic Knowledge Graph"""
    results = []
    if not os.path.exists(kg_path): return None
    import json
    with open(kg_path, 'r', encoding='utf-8') as f:
        kg_data = json.load(f)
        
    if "conditions" not in kg_data: return None
    
    for user_sample in user_data.index:
        sample_expr = user_data.loc[user_sample]
        # Standardize within sample to find the most extreme dysregulation
        z_scores = (sample_expr - sample_expr.mean()) / sample_expr.std()
        # Take the top 250 most dysregulated genes for this sample
        top_genes = set(z_scores.abs().sort_values(ascending=False).head(250).index)
        
        for condition, c_data in kg_data["conditions"].items():
            condition_genes = set([n["id"] for n in c_data["nodes"]])
            overlap = top_genes.intersection(condition_genes)
            activation = len(overlap) / len(condition_genes) if len(condition_genes) > 0 else 0
            
            results.append({
                "User Sample": user_sample,
                "Condition": condition.replace('_', ' '),
                "Activation %": round(activation * 100, 2),
                "Overlap Ratio": f"{len(overlap)} / {len(condition_genes)} nodes",
                "Key Activated Nodes": ", ".join(list(overlap)[:5]) + ("..." if len(overlap) > 5 else "")
            })
            
    return pd.DataFrame(results)

# Main App
st.markdown("### Select Data Source")
col1, col2 = st.columns([2, 1])

with col1:
    uploaded_file = st.file_uploader(
        "Upload Gene Expression Data (CSV or TSV)",
        type=['csv', 'tsv', 'txt'],
        help="Format: Genes as rows, samples as columns. First column should be gene symbols."
    )

with col2:
    st.markdown("<br><br>", unsafe_allow_html=True)
    
    if "use_sample_data" not in st.session_state:
        st.session_state.use_sample_data = False

    if st.button("🚀 Run Sample Data", help="Load sample_expression_data.csv for a quick demonstration"):
        st.session_state.use_sample_data = True

# Also clear the dummy sample state if the user uploads a real file
if uploaded_file:
    st.session_state.use_sample_data = False

if uploaded_file or st.session_state.use_sample_data:
    try:
        # Read file
        if not uploaded_file and st.session_state.use_sample_data:
            user_df = pd.read_csv("sample_expression_data.csv", index_col=0)
            file_name = "sample_expression_data.csv"
            st.success(f"✅ Loaded sample file: {file_name}")
        else:
            if uploaded_file.name.endswith('.csv'):
                user_df = pd.read_csv(uploaded_file, index_col=0)
            else:
                user_df = pd.read_csv(uploaded_file, sep='\t', index_col=0)
            file_name = uploaded_file.name
            st.success(f"✅ Loaded file: {file_name}")
        st.info(f"Dimensions: {user_df.shape[0]} genes × {user_df.shape[1]} samples")
        
        # Show preview
        with st.expander("Preview Data"):
            st.dataframe(user_df.head(10))
        
        # Validation
        validation = validate_user_data(user_df)
        
        if validation['warnings']:
            for warning in validation['warnings']:
                st.warning(f"{warning}")
        
        if not validation['valid']:
            for error in validation['errors']:
                st.error(f"{error}")
            st.stop()
        
        # Preprocessing options
        col1, col2 = st.columns(2)
        with col1:
            auto_log = st.checkbox("Auto log2 transform if raw counts detected", value=True)
        
        # Preprocess
        user_data_processed = preprocess_user_data(user_df, auto_log=auto_log)
        
        st.success(f"✅ Processed data: {user_data_processed.shape[0]} samples × {user_data_processed.shape[1]} genes")
        
        # Load midbase
        with st.spinner("Loading midbase database..."):
            midbase_data, metadata_df = load_midbase_expression()
        
        st.info(f"📚 Midbase: {len(midbase_data)} samples × {len(midbase_data.columns)} genes")
        
        # Calculate similarity
        similarity_results = calculate_similarity(user_data_processed, midbase_data, metadata_df, top_n=20)
        
        if similarity_results is not None:
            # Display results
            st.subheader("Similarity Analysis Results")
            
            # Summary table
            for user_sample in user_data_processed.index:
                st.markdown(f"### Sample: **{user_sample}**")
                
                sample_results = similarity_results[similarity_results['user_sample'] == user_sample].head(5)
                
                # Format table
                display_df = sample_results[['midbase_sample', 'study_id', 'condition', 'correlation']].copy()
                display_df.columns = ['Midbase Sample', 'Study', 'Condition', 'Correlation']
                display_df['Correlation'] = display_df['Correlation'].apply(lambda x: f"{x:.3f}")
                
                st.dataframe(display_df, hide_index=True, width='stretch')
                
                # Top match
                top_match = sample_results.iloc[0]
                st.success(f"**Top Match:** {top_match['condition']} (Study: {top_match['study_id']}, Correlation: {top_match['correlation']:.3f})")
            
            st.markdown("---")
            st.subheader("Interactive Visualizations")
            
            tab1, tab2, tab3, tab4 = st.tabs(["🔵 PCA Overlay", "📊 Full Comparison Data", "🧬 Gene Comparison", "🧠 Graph-Based Prediction"])
            
            with tab1:
                st.markdown("**Principal Component Analysis** - Shows how your data clusters with midbase samples")
                fig_pca = plot_pca_overlay(midbase_data, user_data_processed, metadata_df)
                st.plotly_chart(fig_pca, width='stretch')
                
                with st.expander("📖 **How to Interpret this PCA Graph (Reliability & Guide)**", expanded=True):
                    st.markdown("""
                    **Is PCA a reliable way to compare samples?**
                    Yes, **Principal Component Analysis (PCA)** is a highly reliable statistical method for identifying global patterns and high-level variance across entire transcriptomes. 
                    
                    However, because PCA compresses over 20,000+ genes down into just 2 visual dimensions, it may obscure fine-grained single-gene differences. You should always cross-reference these global clusters with the exact **Pearson Correlation** scores calculated in the summary tables above.
                    
                    **How to read this graph:**
                    - ⭐ **Your Data:** Your uploaded samples are represented by the large **Star** symbols.
                    - 🔵 **Reference Data:** The established midbase diseases are represented by the **Circles**.
                    - **Proximity:** The closer your Star is to a cluster of reference Circles, the more transcriptionally similar your tissue is to that condition.
                    - **Axes (PC1 / PC2):** These axes represent the mathematical directions of highest variance. The percentages (e.g., PC1: 25%) show how much of the total biological complexity is captured by that axis. If your sample sits tightly inside the Klinefelter cluster on both axes, it indicates a very strong global match to that pathology.
                    """)
                
            with tab2:
                st.markdown("**All Comparative Results** - Detailed correlation scores for every sample in the database")
                display_full_df = similarity_results.copy()
                display_full_df.columns = ['Your Sample', 'Midbase Sample', 'Study ID', 'Condition', 'Correlation']
                display_full_df['Correlation'] = display_full_df['Correlation'].apply(lambda x: f"{x:.4f}")
                st.dataframe(display_full_df, hide_index=True, width='stretch')
                
            with tab3:
                st.markdown("**Single Gene Expression Comparison** - Compare your sample against disease averages")
                
                # Allow user to pick genes available in their uploaded dataset
                col_s, col_g = st.columns(2)
                with col_s:
                    selected_sample = st.selectbox("Select Your Sample", options=user_data_processed.index)
                with col_g:
                    # Sort genes alphabetically for easier finding
                    sorted_genes = sorted(user_data_processed.columns.tolist())
                    selected_gene = st.selectbox("Select Gene to Compare", options=sorted_genes)
                
                if selected_sample and selected_gene:
                    user_expr = user_data_processed.loc[selected_sample, selected_gene]
                    
                    # Calculate averages per condition, filtering out Fetal/Stem Cells
                    exclude_conditions = ['Fetal Timeseries', 'Spermatogonial Stem Cells']
                    filtered_meta = metadata_df[~metadata_df['condition'].isin(exclude_conditions)]
                    
                    condition_means = []
                    condition_names = []
                    
                    for condition, group in filtered_meta.groupby('condition'):
                        samples = group['sample_id'].tolist()
                        valid_samples = [s for s in samples if s in midbase_data.index]
                        
                        if valid_samples and selected_gene in midbase_data.columns:
                            mean_expr = midbase_data.loc[valid_samples, selected_gene].mean()
                            condition_means.append(mean_expr)
                            condition_names.append(condition)
                    
                    # Create plotting dataframe
                    plot_df = pd.DataFrame({
                        'Condition': ['Your Sample: ' + str(selected_sample)] + condition_names,
                        'Expression': [user_expr] + condition_means,
                        'Type': ['User Upload'] + ['Reference Database'] * len(condition_names)
                    })
                    
                    # Sort the reference conditions by expression descending, keeping User Upload out first
                    plot_df_ref = plot_df[plot_df['Type'] == 'Reference Database'].sort_values('Expression', ascending=False)
                    plot_df = pd.concat([plot_df[plot_df['Type'] == 'User Upload'], plot_df_ref])
                    
                    fig_bar = px.bar(
                        plot_df,
                        x='Condition',
                        y='Expression',
                        color='Type',
                        color_discrete_map={'User Upload': '#ff2b2b', 'Reference Database': '#005b96'},
                        title=f"{selected_gene} Expression Comparison",
                    )
                    
                    fig_bar.update_layout(xaxis_tickangle=-45)
                    st.plotly_chart(fig_bar, width='stretch')
                    
            with tab4:
                st.markdown("**Graph-Based Disease Prediction** - Evaluates your sample's expression signature against established diagnostic knowledge graphs to predict pathology.")
                
                with st.spinner("Calculating network topology predictions..."):
                    kg_results = calculate_kg_projection(user_data_processed)
                    
                    if kg_results is not None and not len(kg_results) == 0:
                        st.markdown("#### Clinical Predictions")
                        # Generate a top prediction for each user sample
                        for user_sample in kg_results['User Sample'].unique():
                            sample_data = kg_results[kg_results['User Sample'] == user_sample]
                            top_prediction = sample_data.sort_values(by="Activation %", ascending=False).iloc[0]
                            
                            st.success(f"**Sample:** `{user_sample}` ➔ **Predicted Pathology:** **{top_prediction['Condition']}** (Network Activation: {top_prediction['Activation %']}%)")
                        
                        st.markdown("---")
                        st.markdown("**Detailed Topology Overlap Analysis**")
                        # Sort the dataframe so the highest activations are always at the top per sample
                        sorted_kg = kg_results.sort_values(by=["User Sample", "Activation %"], ascending=[True, False])
                        
                        st.dataframe(
                            sorted_kg, 
                            hide_index=True, 
                            use_container_width=True,
                            column_config={
                                "Activation %": st.column_config.ProgressColumn(
                                    "Activation %",
                                    format="%.2f%%",
                                    min_value=0,
                                    max_value=100,
                                )
                            }
                        )
                        st.info("💡 **How this works:** This algorithm identifies the most extremely dysregulated genes in your sample and checks how many of them intersect with the core topological nodes of known diagnostic networks in the Midbase Knowledge Graph.")
                    else:
                        st.warning("⚠️ Could not generate prediction. Ensure `data/knowledge_graph.json` exists securely.")
            
            # Download results
            st.markdown("---")
            st.subheader("Download Results")
            
            csv = similarity_results.to_csv(index=False)
            st.download_button(
                label="Download Similarity Report (CSV)",
                data=csv,
                file_name="similarity_analysis_results.csv",
                mime="text/csv"
            )
    
    except Exception as e:
        st.error(f"Error processing file: {e}")
        import traceback
        st.code(traceback.format_exc())

else:
    # Help section
    st.markdown("---")
    st.subheader("How to Use")
    
    st.markdown("""
    **Step 1:** Prepare your gene expression data
    - Format: CSV or TSV file
    - Genes as rows (first column), samples as columns
    - Gene symbols should match HGNC format (e.g., TP53, BRCA1, AR)
    - Minimum 1,000 genes recommended
    
    **Step 2:** Upload your file
    - System will validate and preprocess automatically
    - Auto-detects if log transformation is needed
    
    **Step 3:** Review similarity results
    - See which disease profiles your data matches
    - Explore interactive visualizations
    - Download detailed report
    
    **Supported Data Types:**
    - Raw counts (will be log2 transformed)
    - Log2-transformed expression
    - TPM/FPKM normalized values
    - DESeq2/edgeR output
    """)
    
    st.markdown("---")
    st.info("**Example:** Upload a CSV with gene expression from testicular tissue to compare against our database of Klinefelter, Azoospermia, AIS, and other conditions.")
