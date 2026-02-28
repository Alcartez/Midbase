"""
Gene Search Page - Query and visualize gene expression across studies
"""
import streamlit as st
import sqlite3
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import requests
import json

st.set_page_config(page_title="Gene Search", page_icon="♂️", layout="wide")

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from style_utils import render_sidebar_nav
render_sidebar_nav("Gene Search")

st.title("Gene Expression Search")

# Database connection
@st.cache_resource
def get_de_db_connection():
    db_path = Path("midbase_de.db")
    return sqlite3.connect(str(db_path), check_same_thread=False)

@st.cache_resource
def get_core_db_connection():
    db_path = Path("midbase_core.db")
    return sqlite3.connect(str(db_path), check_same_thread=False)

@st.cache_data(ttl=3600)
def search_genes(query):
    """Search for genes matching query"""
    conn = get_de_db_connection()
    sql = f"""
    SELECT DISTINCT gene_symbol 
    FROM differential_expression 
    WHERE gene_symbol LIKE '%{query}%'
    ORDER BY gene_symbol
    LIMIT 100
    """
    return pd.read_sql_query(sql, conn)

@st.cache_data(ttl=3600)
def get_gene_de_results(gene_symbol):
    """Get all DE results for a specific gene by manually joining core and DE databases"""
    conn_de = get_de_db_connection()
    conn_core = get_core_db_connection()
    
    # 1. Fetch DE data
    sql_de = f"""
    SELECT 
        study_id,
        gene_symbol,
        logFC,
        ave_expr,
        p_value,
        comparison
    FROM differential_expression 
    WHERE gene_symbol = '{gene_symbol}'
    """
    df_de = pd.read_sql_query(sql_de, conn_de)
    
    if df_de.empty:
        return df_de
        
    # 2. Fetch Studies metadata
    study_ids = tuple(df_de['study_id'].unique())
    if len(study_ids) == 1:
        # SQLite syntax for single tuple IN clause requires trailing comma
        study_ids_str = f"('{study_ids[0]}')"
    else:
        study_ids_str = str(study_ids)
        
    sql_studies = f"SELECT study_id, title as study_title FROM studies WHERE study_id IN {study_ids_str}"
    df_studies = pd.read_sql_query(sql_studies, conn_core)
    
    # 3. Merge in pandas
    df_merged = pd.merge(df_de, df_studies, on='study_id', how='left')
    df_merged = df_merged.sort_values('p_value')
    
    return df_merged

@st.cache_resource
def get_meta_db_connection():
    db_path = Path("gene_metadata.db")
    if db_path.exists():
        return sqlite3.connect(str(db_path), check_same_thread=False)
    return None

@st.cache_data(ttl=3600)
def fetch_gene_info(gene_symbol):
    """Fetch gene details from local SQLite metadata cache first, fallback to API if needed"""
    conn = get_meta_db_connection()
    if conn:
        try:
            sql = f"SELECT name, type, aliases, summary, kegg_pathways, go_bp FROM gene_info WHERE gene_symbol = ?"
            df = pd.read_sql_query(sql, conn, params=(gene_symbol,))
            
            if not df.empty:
                row = df.iloc[0]
                
                # Parse JSON lists back to Python lists
                import json
                aliases = json.loads(row['aliases']) if row['aliases'] else []
                kegg = json.loads(row['kegg_pathways']) if row['kegg_pathways'] else []
                go_bp = json.loads(row['go_bp']) if row['go_bp'] else []
                
                return {
                    'name': row['name'] or '',
                    'summary': row['summary'] or '',
                    'aliases': aliases,
                    'type': row['type'] or 'Unknown',
                    'kegg': kegg[:5],  # Keep reasonable limits for UI display
                    'go_bp': go_bp[:6]
                }
        except Exception as e:
            st.warning(f"Error reading local metadata cache: {e}")
            
    # Fallback to remote API if DB missing or gene not found
    try:
        url = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=name,summary,alias,type_of_gene,MIM,pathway.kegg,go.BP&size=1"
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            data = r.json()
            if data.get('hits'):
                hit = data['hits'][0]
                
                # Basic parsing
                aliases = hit.get('alias', [])
                if isinstance(aliases, str): aliases = [aliases]
                
                kegg_paths = []
                pathways = hit.get('pathway', {}).get('kegg', []) if isinstance(hit.get('pathway'), dict) else []
                if isinstance(pathways, dict): pathways = [pathways]
                for p in (pathways or [])[:5]:
                    if isinstance(p, dict) and 'name' in p: kegg_paths.append(p['name'])
                
                go_terms = []
                go_pb = hit.get('go', {}).get('BP', []) if isinstance(hit.get('go'), dict) else []
                if isinstance(go_pb, dict): go_pb = [go_pb]
                for g in (go_pb or [])[:6]:
                    if isinstance(g, dict) and 'term' in g: go_terms.append(g['term'])
                
                return {
                    'name': hit.get('name', ''),
                    'summary': hit.get('summary', ''),
                    'aliases': aliases,
                    'type': hit.get('type_of_gene', 'Unknown'),
                    'kegg': kegg_paths,
                    'go_bp': go_terms,
                }
    except Exception:
        pass
        
    return None

def build_gene_context(gene_symbol, gene_info, de_results):
    """Build a rich context string for the LLM"""
    ctx = f"=== Gene: {gene_symbol} ===\n"
    if gene_info:
        ctx += f"Full name: {gene_info.get('name','')}\n"
        if gene_info.get('summary'):
            ctx += f"Summary: {gene_info['summary'][:600]}\n"
        if gene_info.get('kegg'):
            ctx += f"KEGG pathways: {', '.join(gene_info['kegg'])}\n"
        if gene_info.get('go_bp'):
            ctx += f"Biological processes: {', '.join(gene_info['go_bp'])}\n"
    ctx += f"\n=== Differential Expression in Testicular Studies ===\n"
    ctx += de_results[['study_id','logFC','p_value','comparison']].to_string(index=False)
    return ctx

# Provider configurations
PROVIDERS = {
    "Groq (Free)": {
        "url": "https://api.groq.com/openai/v1",
        "model": "llama-3.3-70b-versatile",
        "mode": "openai",
        "help": "Free tier — get key at console.groq.com",
    },
    "Google Gemini": {
        "url": "https://generativelanguage.googleapis.com",
        "model": "gemini-1.5-flash",
        "mode": "gemini",
        "help": "Free tier — get key at aistudio.google.com",
    },
    "OpenAI": {
        "url": "https://api.openai.com/v1",
        "model": "gpt-4o-mini",
        "mode": "openai",
        "help": "Paid — get key at platform.openai.com",
    },
    "Anthropic Claude": {
        "url": "https://api.anthropic.com/v1",
        "model": "claude-3-haiku-20240307",
        "mode": "claude",
        "help": "Paid — get key at console.anthropic.com",
    },
}

def chat_with_llm(messages, provider_cfg, api_key):
    """Dispatch to correct LLM SDK/API based on provider mode"""
    mode  = provider_cfg["mode"]
    model = provider_cfg["model"]

    if mode == "gemini":
        try:
            import google.generativeai as genai
            genai.configure(api_key=api_key)
            # Convert messages: skip system role (unsupported in basic API)
            history = []
            system_txt = ""
            for m in messages:
                if m["role"] == "system":
                    system_txt = m["content"]
                elif m["role"] == "user":
                    history.append(m["content"])
            full_prompt = (system_txt + "\n\n" if system_txt else "") + (history[-1] if history else "")
            g_model = genai.GenerativeModel(model)
            resp = g_model.generate_content(full_prompt)
            return resp.text
        except ImportError:
            raise RuntimeError("google-generativeai not installed. Run: pip install google-generativeai")

    elif mode == "claude":
        headers = {
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        }
        # Extract system message
        system_msg = next((m["content"] for m in messages if m["role"] == "system"), "")
        non_system = [m for m in messages if m["role"] != "system"]
        payload = {
            "model": model,
            "max_tokens": 512,
            "system": system_msg,
            "messages": non_system,
        }
        r = requests.post("https://api.anthropic.com/v1/messages",
                          headers=headers, json=payload, timeout=30)
        r.raise_for_status()
        return r.json()["content"][0]["text"]

    else:  # openai-compatible (Groq, OpenAI, Ollama, OpenRouter)
        headers = {"Content-Type": "application/json",
                   "Authorization": f"Bearer {api_key}"}
        payload = {"model": model, "messages": messages,
                   "max_tokens": 512, "temperature": 0.7}
        r = requests.post(f"{provider_cfg['url'].rstrip('/')}/chat/completions",
                          headers=headers, json=payload, timeout=30)
        r.raise_for_status()
        return r.json()["choices"][0]["message"]["content"]

# Setup state
current_gene = st.session_state.get('selected_gene')
gene_info = None

if current_gene:
    with st.spinner(f"Fetching {current_gene} info..."):
        gene_info = fetch_gene_info(current_gene)

# ---- AI Assistant Sidebar ----
with st.sidebar:
    st.header("Chat Context")

    with st.expander("Provider Settings", expanded=False):
        provider_name = st.selectbox("Provider", list(PROVIDERS.keys()))
        provider_cfg  = PROVIDERS[provider_name].copy()

        st.caption(provider_cfg["help"])

        api_key = st.text_input("API Key", type="password", key="llm_api_key")

        # Allow overriding the default model
        provider_cfg["model"] = st.text_input(
            "Model", value=provider_cfg["model"], key="llm_model"
        )
        llm_ready = bool(api_key)

    st.markdown("---")

    if current_gene:
        st.markdown(f"**Chatting about: {current_gene}**")

        # Chat
        if 'chat_messages' not in st.session_state or st.session_state.get('chat_gene') != current_gene:
            # Get DE results for system context
            de_for_ctx = get_gene_de_results(current_gene)
            gene_ctx = build_gene_context(current_gene, gene_info if current_gene else None, de_for_ctx)
            st.session_state.chat_messages = [
                {"role": "system", "content":
                    f"You are a testicular biology expert assistant. Answer concisely (2-4 sentences). "
                    f"Use the following context about the current gene:\n\n{gene_ctx}"}
            ]
            st.session_state.chat_gene = current_gene
            st.session_state.chat_history = []  # displayed history

        # Show chat history
        for msg in st.session_state.get('chat_history', []):
            role_icon = "" if msg['role'] == 'user' else ""
            st.markdown(f"{role_icon} {msg['content']}")

        user_q = st.text_area("", placeholder=f"e.g. Is {current_gene} dysregulated in azoospermia?",
                              height=70, label_visibility='collapsed', key="chat_input")

        col_ask, col_clear = st.columns([2, 1])
        with col_ask:
            ask_btn = st.button("Ask", use_container_width=True, type="primary")
        with col_clear:
            if st.button("Clear", use_container_width=True):
                st.session_state.chat_history = []
                st.session_state.chat_messages = st.session_state.chat_messages[:1]
                st.rerun()

        if ask_btn and user_q.strip():
            if not llm_ready:
                st.warning("Enter your API key above to enable chat")
            else:
                st.session_state.chat_messages.append({"role": "user", "content": user_q})
                st.session_state.chat_history.append({"role": "user", "content": user_q})
                with st.spinner(f"Asking {provider_name}..."):
                    try:
                        reply = chat_with_llm(
                            st.session_state.chat_messages,
                            provider_cfg, api_key
                        )
                        st.session_state.chat_messages.append({"role": "assistant", "content": reply})
                        st.session_state.chat_history.append({"role": "assistant", "content": reply})
                    except Exception as e:
                        st.error(f"LLM error: {e}")
                st.rerun()
    else:
        st.info("Search for a gene to see its details and chat about it")

# Search interface
col_search, col_filter = st.columns([3, 1])

# Initialize session state for selected gene
if 'selected_gene' not in st.session_state:
    st.session_state.selected_gene = None

with col_search:
    gene_input = st.text_input(
        "Enter gene symbol",
        placeholder="e.g., TP53, BRCA1, AR, SOX9",
        help="Search for gene expression and differential expression results"
    )

with col_filter:
    st.markdown("<br>", unsafe_allow_html=True)
    search_btn = st.button("Search", use_container_width=True, type="primary")

# Trigger search on button click
if search_btn and gene_input:
    st.session_state.selected_gene = gene_input.upper()

# Display results if a gene is selected
if st.session_state.selected_gene:
    gene_to_search = st.session_state.selected_gene
    
    # Get DE results
    de_results = get_gene_de_results(gene_to_search)
    
    if len(de_results) > 0:
        st.success(f"Found gene: **{gene_to_search}**")
        
        # Display Gene Info Card in Main Area
        if gene_info:
            with st.container():
                st.markdown(f"""
                <div style='background-color: white; border: 1px solid #e1e8ed; border-radius: 8px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.02);'>
                    <h3 style='margin-top: 0; color: #005b96;'>{gene_to_search} — {gene_info.get('name', '')}</h3>
                    <div style='color: #64748b; font-size: 0.9em; margin-bottom: 12px;'>
                        <b>Type:</b> {gene_info.get('type', 'Unknown')} &nbsp;|&nbsp; 
                        <b>Aliases:</b> {', '.join(gene_info.get('aliases', [])[:6]) if gene_info.get('aliases') else 'None'}
                    </div>
                """, unsafe_allow_html=True)
                
                if gene_info.get('summary'):
                    st.markdown(f"**Summary:** {gene_info['summary']}")
                
                col_kegg, col_go = st.columns(2)
                with col_kegg:
                    if gene_info.get('kegg'):
                        st.markdown("**Pathways (KEGG):**")
                        for p in gene_info['kegg']:
                            st.markdown(f"- {p}")
                with col_go:
                    if gene_info.get('go_bp'):
                        st.markdown("**Biological Processes (GO):**")
                        for p in gene_info['go_bp']:
                            st.markdown(f"- {p}")
                
                st.markdown("</div>", unsafe_allow_html=True)
        else:
            st.info("Gene details (functions/pathways) unavailable from MyGene.info at this moment.")
            
        st.markdown("---")
        
        # Tabs for different views
        tab1, tab2, tab3 = st.tabs(["Summary", "Volcano Plot", "Data Table"])
        
        with tab1:
            st.subheader(f"Differential Expression Summary - {gene_to_search}")
            
            # Stats cards
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Studies", de_results['study_id'].nunique())
            with col2:
                sig_count = (de_results['p_value'] < 0.05).sum()
                st.metric("Significant (p<0.05)", sig_count)
            with col3:
                max_logfc = de_results['logFC'].abs().max()
                st.metric("Max |logFC|", f"{max_logfc:.2f}")
            
            st.markdown("")
            
            # Bar chart of logFC across studies
            fig = px.bar(
                de_results,
                x='study_id',
                y='logFC',
                color='p_value',
                color_continuous_scale='RdBu_r',
                title=f"{gene_to_search} - Log Fold Change Across Studies",
                labels={'logFC': 'Log2 Fold Change', 'study_id': 'Study', 'p_value': 'P-value'},
                hover_data=['study_title', 'comparison', 'ave_expr']
            )
            fig.add_hline(y=0, line_dash="dash", line_color="gray")
            fig.update_layout(height=400)
            st.plotly_chart(fig, width='stretch')
        
        with tab2:
            st.subheader("Volcano Plot Visualization")
            
            # Use p_value directly
            
            # Volcano plot (using p-value directly)
            fig = px.scatter(
                de_results,
                x='logFC',
                y='p_value',
                hover_name='study_id',
                hover_data=['comparison'],
                labels={'logFC': 'Log2 Fold Change', 'p_value': 'P-value'},
                title=f'{gene_to_search} - Differential Expression Across Studies'
            )
            
            # Significance threshold
            fig.add_hline(y=0.05, line_dash="dash", line_color="red", 
                         annotation_text="p=0.05")
            fig.add_vline(x=1, line_dash="dash", line_color="gray")
            fig.add_vline(x=-1, line_dash="dash", line_color="gray")
            
            fig.update_layout(height=500)
            st.plotly_chart(fig, width='stretch')
        
        with tab3:
            st.subheader("Detailed Results Table")
            
            # Format p-values
            de_results_display = de_results.copy()
            de_results_display['p_value'] = de_results_display['p_value'].apply(lambda x: f"{x:.2e}")
            de_results_display['logFC'] = de_results_display['logFC'].apply(lambda x: f"{x:.3f}")
            de_results_display['ave_expr'] = de_results_display['ave_expr'].apply(lambda x: f"{x:.2f}")
            
            st.dataframe(
                de_results_display,
                width='stretch',
                column_config={
                    "study_id": "Study ID",
                    "study_title": "Study Title",
                    "gene_symbol": "Gene",
                    "logFC": "Log Fold Change",
                    "ave_expr": "Avg Expression",
                    "p_value": "P-value",
                    "comparison": "Comparison"
                }
            )
            
            # Download button
            csv = de_results.to_csv(index=False)
            st.download_button(
                label="Download Results (CSV)",
                data=csv,
                file_name=f"{gene_to_search}_DE_results.csv",
                mime="text/csv"
            )
        
        # Clear button
        if st.button("🔄 Search Another Gene"):
            st.session_state.selected_gene = None
            st.rerun()
    else:
        st.warning(f"Gene '{gene_to_search}' not found in database.")
        if st.button("Try Another Gene"):
            st.session_state.selected_gene = None
            st.rerun()

elif search_btn and gene_input:
    # Search for matching genes only if button clicked and no exact match
    matches = search_genes(gene_input.upper())
    
    if len(matches) == 0:
        st.warning(f"No genes found matching '{gene_input}'")
    elif len(matches) == 1:
        # Exact match - set it
        st.session_state.selected_gene = matches.iloc[0]['gene_symbol']
        st.rerun()
    else:
        # Multiple matches - show selection
        st.info(f"Found {len(matches)} genes matching '{gene_input}'. Select one:")
        
        # Create selection columns
        cols = st.columns(5)
        for idx, row in matches.iterrows():
            col_idx = idx % 5
            with cols[col_idx]:
                if st.button(row['gene_symbol'], key=f"gene_{idx}"):
                    st.session_state.selected_gene = row['gene_symbol']
                    st.rerun()

# Help section
if not gene_input or not search_btn:
    st.markdown("---")
    st.subheader("💡 How to Use")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        **Search Tips:**
        - Enter complete gene symbols (e.g., TP53, BRCA1)
        - Partial matches supported (e.g., "SOX" finds SOX9, SOX2, etc.)
        - Case insensitive
        
        **Interpreting Results:**
        - **logFC > 0**: Gene is upregulated in disease
        - **logFC < 0**: Gene is downregulated in disease
        - **p < 0.05**: Statistically significant
        """)
    
    with col2:
        st.markdown("""
        **Example Genes to Try:**
        - **TP53**: Tumor suppressor
        - **AR**: Androgen receptor
        - **SOX9**: Sex determination
        - **DMRT1**: Male development
        - **AMH**: Anti-Müllerian hormone
        """)
