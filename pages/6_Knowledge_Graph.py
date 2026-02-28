import streamlit as st
import json
import time
import pandas as pd
import numpy as np
from pathlib import Path
import sys
from streamlit_agraph import agraph, Node, Edge, Config
sys.path.append(str(Path(__file__).parent.parent))
from style_utils import render_sidebar_nav

st.set_page_config(
    page_title="Diagnostic Knowledge Graph",
    page_icon="♂️",
    layout="wide"
)

render_sidebar_nav("Knowledge Graph")

st.title("Diagnostic Knowledge Graph Explorer")

st.markdown("""
This interactive network visualizes the core expression and physical protein interactions defining major male infertility conditions. 
* **Green Edges:** Known physical Protein-Protein Interactions (STRING DB).
* **Gray Edges:** Strong co-expression observed within our clinical samples (|R| > threshold).
""")


# --- Data Loading with Custom UI ---

@st.cache_data(show_spinner=False, ttl=3600)
def fetch_knowledge_graph_json():
    """Load the precomputed knowledge graph strictly from disk with no caching visible spinner."""
    # Simulate a tiny load time if the user is perceiving the 'proper loading screen' (user request)
    time.sleep(1.5) 
    kr_path = Path("data/knowledge_graph.json")
    if not kr_path.exists():
        return None
    with open(kr_path, "r", encoding="utf-8") as f:
        return json.load(f)

# The user explicitly asked for a proper loading screen, NOT the default Streamlit spinner.
# We will inject a sleek custom HTML/CSS CSS loader container, then remove it once loaded.
if "kg_loaded" not in st.session_state:
    st.session_state.kg_loaded = False

if not st.session_state.kg_loaded:
    loader_placeholder = st.empty()
    with loader_placeholder.container():
        st.markdown("""
            <style>
            .custom-loader-wrapper {
                display: flex;
                flex-direction: column;
                align-items: center;
                justify-content: center;
                height: 400px;
                background: rgba(255, 255, 255, 0.05);
                border-radius: 15px;
                border: 1px solid rgba(255, 255, 255, 0.1);
            }
            .spinner {
                width: 50px;
                height: 50px;
                border: 5px solid rgba(0, 150, 255, 0.2);
                border-top: 5px solid #0096FF;
                border-radius: 50%;
                animation: spin 1s linear infinite;
                margin-bottom: 20px;
            }
            @keyframes spin {
                0% { transform: rotate(0deg); }
                100% { transform: rotate(360deg); }
            }
            .loading-text {
                font-family: 'Inter', sans-serif;
                font-size: 1.2rem;
                color: #A0AEC0;
                font-weight: 500;
            }
            </style>
            <div class="custom-loader-wrapper">
                <div class="spinner"></div>
                <div class="loading-text">Assembling Multi-Dimensional Subgraphs...</div>
            </div>
            """, unsafe_allow_html=True)
        
    # Trigger actual data load while the custom UI is showing
    kg_data = fetch_knowledge_graph_json()
    st.session_state.kg_data = kg_data
    st.session_state.kg_loaded = True
    
    # Rerun to clear the loader and show the graph
    st.rerun()

kg_data = st.session_state.kg_data

if not kg_data:
    st.error("Knowledge Graph data not found. Please run scripts/build_knowledge_graph.py first.")
    st.stop()


# --- Main App Logic ---

conditions = list(kg_data["conditions"].keys())

col1, col2 = st.columns([1, 3])

with col1:
    st.subheader("Network Filters")
    selected_condition = st.selectbox("Select Diagnostic Profile", conditions)
    
    # Filter edges by co-expression strictness
    st.markdown("---")
    co_expr_thresh = st.slider("Co-expression Threshold (|R|)", min_value=0.70, max_value=0.99, value=0.85, step=0.01,
                              help="Increase to show only the strongest dataset-specific co-expression links (gray edges).")
    max_edges = st.slider("Max Visible Edges", min_value=50, max_value=1200, value=350, step=50,
                              help="Limit the number of edges. If physics glitch, reduce this number.")
    spring_length = st.slider("Edge Spring Length", min_value=50, max_value=350, value=120, step=10,
                              help="Increase this value to push nodes further apart and stop 'dancing graphs'.")
    
    st.markdown("---")
    st.markdown("### Edge Legend")
    st.markdown("🟢 **Green:** Physical Protein Interaction (PPI)")
    st.markdown("🔵 **Blue:** Functional Gene Interaction (GGI)")
    st.markdown("⚪ **Gray:** Clinical Co-expression")
    
# Process the graph for the selected condition
with col2:
    condition_data = kg_data["conditions"][selected_condition]
    
    nodes_data = condition_data["nodes"]
    edges_data = condition_data["edges"]
    
    nodes = []
    edges = []
    
    # We need to map nodes clearly. For agraph, we append Node()
    for n in nodes_data:
        # Color based on Fold Change direction (LogFC)
        # Red = Upregulated, Blue = Downregulated (standard DE colors)
        logfc = n.get("logFC", 0)
        color = "#ff4b4b" if logfc > 0 else "#2b83ff"
        
        nodes.append(Node(
            id=n["id"],
            label=n["label"],
            size=25,
            color=color,
            title=f"Gene: {n['label']}\nLogFC: {logfc}\np-value: {n.get('p_value', 0):.2e}"
        ))
        
    # Pre-filter edges to prevent browser freezing (agraph struggles > 500)
    valid_edges = []
    for e in edges_data:
        weight = float(e.get("weight", 0))
        if e["type"] == "coexpression" and abs(weight) < co_expr_thresh:
            continue
        valid_edges.append(e)
        
    # Cap total edges by keeping only the strongest connections
    # Prioritize biological STRING edges (PPI/GGI) so they aren't hidden by hundreds of co-expression edges
    def edge_sort_key(x):
        type_bonus = 10.0 if x.get("type") in ["ppi", "ggi"] else 0.0
        return type_bonus + abs(float(x.get("weight", 0)))
        
    if len(valid_edges) > max_edges:
        valid_edges = sorted(valid_edges, key=edge_sort_key, reverse=True)[:max_edges]
        
    for e in valid_edges:
        weight = float(e.get("weight", 0))
        edges.append(Edge(
            source=e["source"],
            target=e["target"],
            color=e["color"],
            title=f"{e['type'].upper()} Weight: {weight:.2f}",
            width=2 if e["type"] in ["ppi", "ggi"] else 1
        ))

    # Configure graph layout and physics
    config = Config(
        width='100%',
        height=600,
        directed=False, 
        physics=True,
        hierarchical=False,
        nodeHighlightBehavior=True,
        highlightColor="#F7A7A6",
        collapsible=False,
        # Enable physics stabilization so it doesn't look completely crazy on load
        # We can pass kwargs to access full vis.js physics config
        physics_kwargs={
            "barnesHut": {
                "gravitationalConstant": -40000,
                "centralGravity": 0.3,
                "springLength": spring_length,
                "springConstant": 0.04,
                "damping": 0.09,
                "avoidOverlap": 1.0
            }
        }
    )

import plotly.express as px
import sqlite3
import requests

@st.cache_data(show_spinner=False, ttl=3600)
def fetch_network_pathways(gene_list):
    """Run Enrichr on the network's genes via REST API to bypass scipy DLL issues"""
    if not gene_list: return None
    try:
        # 1. Add list to Enrichr
        add_url = "https://maayanlab.cloud/Enrichr/addList"
        payload = {
            'list': (None, '\n'.join(gene_list)),
            'description': (None, 'KG_Network_Genes')
        }
        res = requests.post(add_url, files=payload)
        if not res.ok: return None
        data = res.json()
        user_list_id = data['userListId']
        
        # 2. Get enrichment results
        libraries = ['KEGG_2021_Human', 'GO_Biological_Process_2023']
        all_results = []
        
        for lib in libraries:
            enrich_url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={lib}"
            res_enrich = requests.get(enrich_url)
            if not res_enrich.ok: continue
            
            enrich_data = res_enrich.json()
            if lib not in enrich_data: continue
            
            for row in enrich_data[lib]:
                term = row[1]
                adj_p = row[6]
                if adj_p < 0.05:
                    all_results.append({
                        'Term': term,
                        'Adjusted P-value': adj_p,
                        'Gene_set': lib
                    })
                    
        if not all_results: return None
        
        results_df = pd.DataFrame(all_results)
        return results_df.sort_values('Adjusted P-value').head(12)
        
    except Exception as e:
        print(f"Enrichr API Error: {e}")
        return None

def fetch_local_gene_info(gene):
    """Fetch gene details from our local precomputed DB"""
    try:
        conn = sqlite3.connect("gene_metadata.db")
        sql = "SELECT name, summary, aliases FROM gene_info WHERE gene_symbol = ?"
        df = pd.read_sql_query(sql, conn, params=(gene,))
        conn.close()
        if not df.empty:
            row = df.iloc[0]
            aliases = json.loads(row['aliases']) if row['aliases'] else []
            return {"name": row['name'], "summary": row['summary'], "aliases": aliases}
    except:
        pass
    return None

# Calculate pathways while graph renders
graph_genes = [n.id for n in nodes]
pathway_df = fetch_network_pathways(graph_genes)

with col2:
    st.subheader(f"Subgraph Mapping: {selected_condition.replace('_', ' ')}")
    
    if len(nodes) > 0:
        actual_edges = len(edges)
        st.caption(f"Visualizing {len(nodes)} core genes with {actual_edges} interactions.")
    
    with st.container(border=True):
        clicked_node = agraph(nodes=nodes, edges=edges, config=config)

    # Display selection details instantly if a node is clicked
    if clicked_node:
        with st.expander(f"🧬 Gene Focus: {clicked_node}", expanded=True):
            info = fetch_local_gene_info(clicked_node)
            if info and info.get('name'):
                st.markdown(f"**Full Name:** {info['name']}")
                if info.get('aliases'):
                    st.caption(f"Aliases: {', '.join(info['aliases'][:5])}")
                st.markdown(info['summary'])
            else:
                st.info("No extended summary found in local cache.")

# --- Pathway Enrichment Section ---
st.markdown("---")
st.subheader("Network External Pathways & Ontology")

if pathway_df is not None and not pathway_df.empty:
    # Format for plotting
    pathway_df['Term'] = pathway_df['Term'].apply(lambda x: (x[:45] + '...') if len(x) > 45 else x)
    pathway_df['-log10(P-value)'] = -np.log10(pathway_df['Adjusted P-value'])
    
    fig = px.bar(
        pathway_df.sort_values('-log10(P-value)', ascending=True),
        x='-log10(P-value)', 
        y='Term',
        color='Gene_set',
        orientation='h',
        title=f"Top Enriched Pathways for the {selected_condition.replace('_', ' ')} Network",
        labels={'-log10(P-value)': 'Significance (-log10 P.Adj)'},
        color_discrete_map={"KEGG_2021_Human": "#4C78A8", "GO_Biological_Process_2023": "#F58518"}
    )
    
    fig.update_layout(height=400, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("No significantly enriched pathways found for this gene network.")

