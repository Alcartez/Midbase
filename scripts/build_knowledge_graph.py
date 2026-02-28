import os
import sqlite3
import pandas as pd
import numpy as np
import requests
import json
from pathlib import Path

# Paths
# Paths
DB_PATH = Path("../midbase_de.db")
PROCESSED_DATA_DIR = Path("../data/processed")
OUTPUT_KG_PATH = Path("../data/knowledge_graph.json")

# STRING API Settings
STRING_API_URL = "https://string-db.org/api"
OUTPUT_FORMAT = "json"
METHOD = "network"

# Conditions to profile (must match study_id or condition mappings)
# Let's target the major ones we have DE results for
TARGET_CONDITIONS = {
    "GSE103905": "Klinefelter_Syndrome",
    "GSE106487_OA": "Obstructive_Azoospermia",
    "GSE106487_NOA": "Non_Obstructive_Azoospermia",
    "GSE125222": "Androgen_Insensitivity_Syndrome",
    "GSE235210": "Zika_Virus_Infection"
}

def get_top_de_genes(study_id, limit=50):
    """Fetch top highly significant DE genes for a study/condition"""
    print(f"Fetching top {limit} DE genes for {study_id}...")
    conn = sqlite3.connect(DB_PATH)
    
    # Handle the cross-study names if necessary
    if study_id in ['GSE106487_OA', 'GSE106487_NOA']:
        base_study = "GSE106487"
        comp_target = "Obstructive_Azoospermia" if "OA" in study_id and "NOA" not in study_id else "Non_Obstructive_Azoospermia"
        sql = f"""
            SELECT gene_symbol, logFC, p_value 
            FROM differential_expression 
            WHERE study_id = '{base_study}' AND comparison LIKE '%vs_{comp_target}%' AND p_value < 0.05
            ORDER BY ABS(logFC) DESC
            LIMIT {limit}
        """
    else:
        sql = f"""
            SELECT gene_symbol, logFC, p_value 
            FROM differential_expression 
            WHERE study_id = '{study_id}' AND p_value < 0.05
            ORDER BY ABS(logFC) DESC
            LIMIT {limit}
        """
        
    df = pd.read_sql_query(sql, conn)
    conn.close()
    return df

def compute_coexpression_edges(study_id, gene_list, threshold=0.8):
    """Read the raw expression matrix and compute gene-gene Pearson correlation"""
    print(f"Computing co-expression edges for {len(gene_list)} genes...")
    
    # Determine which file to load
    base_study = study_id.split('_')[0]
    csv_path = PROCESSED_DATA_DIR / f"{base_study}.csv"
    
    if not csv_path.exists():
        print(f"  Warning: Raw data {csv_path} not found. Skipping co-expression.")
        return []
        
    df = pd.read_csv(csv_path, index_col=0)
    
    # Optional: Filter columns to only the condition samples if it's a mixed dataset
    # e.g., for GSE106487_OA, only use OA columns. For simplicity, we can use the whole matrix 
    # if it represents the diseased state relative to normal, but disease-only is better.
    if study_id == "GSE106487_OA":
        cols = [c for c in df.columns if 'Obstructive' in c and 'Non' not in c]
        if cols: df = df[cols]
    elif study_id == "GSE106487_NOA":
        cols = [c for c in df.columns if 'Non_Obstructive' in c]
        if cols: df = df[cols]
    
    # Filter matrix to only our target genes
    available_genes = [g for g in gene_list if g in df.columns]
    expr_matrix = df[available_genes] # Samples are rows, genes are columns
    
    # Calculate Pearson Correlation
    corr_matrix = expr_matrix.corr(method='pearson')
    
    edges = []
    # Extract upper triangle (avoid duplicates and self-loops)
    for i in range(len(corr_matrix.columns)):
        for j in range(i+1, len(corr_matrix.columns)):
            gene_a = corr_matrix.columns[i]
            gene_b = corr_matrix.columns[j]
            weight = corr_matrix.iloc[i, j]
            
            if abs(weight) >= threshold:
                edges.append({
                    "source": gene_a,
                    "target": gene_b,
                    "type": "coexpression",
                    "weight": round(weight, 3),
                    "color": "gray"
                })
                
    print(f"  Found {len(edges)} co-expression edges (|R| >= {threshold})")
    return edges

def fetch_string_edges(gene_list, required_score=400):
    """Fetch physical (PPI) and functional (GGI) interactions from STRING DB API"""
    print(f"Fetching STRING physical and functional edges for {len(gene_list)} genes...")
    
    request_url = f"{STRING_API_URL}/{OUTPUT_FORMAT}/{METHOD}"
    ppi_params = {
        "identifiers": "%0d".join(gene_list),
        "species": 9606,
        "required_score": required_score,
        "network_type": "physical",
        "caller_identity": "xy_counsel_app"
    }
    
    ggi_params = {
        "identifiers": "%0d".join(gene_list),
        "species": 9606,
        "required_score": required_score,
        "network_type": "functional",
        "caller_identity": "xy_counsel_app"
    }
    
    edges = []
    physical_pairs = set()
    
    try:
        # 1. Fetch Physical PPIs
        res_phys = requests.post(request_url, data=ppi_params, timeout=15)
        res_phys.raise_for_status()
        
        for row in res_phys.json():
            source = row["preferredName_A"]
            target = row["preferredName_B"]
            pair = tuple(sorted([source, target]))
            physical_pairs.add(pair)
            edges.append({
                "source": source,
                "target": target,
                "type": "ppi",
                "weight": row["score"],
                "color": "green"
            })
            
        print(f"  Found {len(physical_pairs)} STRING PPI (Physical) edges")
        
        # 2. Fetch Functional GGIs
        res_func = requests.post(request_url, data=ggi_params, timeout=15)
        res_func.raise_for_status()
        
        ggi_count = 0
        for row in res_func.json():
            source = row["preferredName_A"]
            target = row["preferredName_B"]
            pair = tuple(sorted([source, target]))
            
            # If it is a functional link but NOT physical, it's a pure GGI
            if pair not in physical_pairs:
                edges.append({
                    "source": source,
                    "target": target,
                    "type": "ggi",
                    "weight": row["score"],
                    "color": "#3b82f6"  # distinct blue color for Gene-Gene Interaction
                })
                physical_pairs.add(pair)
                ggi_count += 1
                
        print(f"  Found {ggi_count} STRING GGI (Functional) edges")
        return edges
        
    except Exception as e:
        print(f"  Error fetching STRING data: {e}")
        return []

def build_knowledge_graph():
    """Main pipeline to construct the multi-layered Knowledge Graph"""
    print("=== Building Diagnostic Knowledge Graph ===")
    
    kg_data = {
        "metadata": {
            "description": "XY Counsel Diagnostic Knowledge Graph",
            "edge_types": {
                "coexpression": {"color": "gray", "description": "High Pearson Correlation (|R| > 0.8) within patient class"},
                "ppi": {"color": "green", "description": "STRING Protein-Protein Interaction (Physical)"},
                "ggi": {"color": "#3b82f6", "description": "STRING Gene-Gene Interaction (Functional)"}
            }
        },
        "conditions": {}
    }
    
    for study_id, condition_name in TARGET_CONDITIONS.items():
        print(f"\\nProcessing: {condition_name}")
        
        # 1. Get nodes
        de_genes = get_top_de_genes(study_id, limit=60) # Top 60 for visual density
        
        if de_genes.empty:
            print("  No DE genes found. Skipping.")
            continue
            
        gene_list = de_genes['gene_symbol'].tolist()
        
        # Build node dictionaries
        nodes = []
        for _, row in de_genes.iterrows():
            nodes.append({
                "id": row['gene_symbol'],
                "label": row['gene_symbol'],
                "logFC": round(row['logFC'], 2),
                "p_value": float(row['p_value']) # JSON serializeable
            })
            
        # 2. Compute Internal Co-expression Edges
        coexpr_edges = compute_coexpression_edges(study_id, gene_list, threshold=0.7) # Allow slightly more leniency for sparse datasets
        
        # 3. Fetch External PPI and GGI Edges
        string_edges = fetch_string_edges(gene_list, required_score=400) # Medium confidence
        
        # 4. Save to master JSON
        kg_data["conditions"][condition_name] = {
            "nodes": nodes,
            "edges": coexpr_edges + string_edges
        }
        
    # Write to disk
    os.makedirs(OUTPUT_KG_PATH.parent, exist_ok=True)
    with open(OUTPUT_KG_PATH, 'w', encoding='utf-8') as f:
        json.dump(kg_data, f, indent=2)
        
    print(f"\\n✅ Knowledge Graph saved to {OUTPUT_KG_PATH}")
    print(f"Total conditions mapped: {len(kg_data['conditions'])}")

if __name__ == "__main__":
    # Ensure run from right dir if relative paths are used
    if not os.path.exists("../midbase_de.db") and os.path.exists("midbase_de.db"):
        DB_PATH = Path("midbase_de.db")
        PROCESSED_DATA_DIR = Path("data/processed")
        OUTPUT_KG_PATH = Path("data/knowledge_graph.json")
        
    build_knowledge_graph()
