import sqlite3
import pandas as pd
import requests
import json
import time

def setup_db():
    conn = sqlite3.connect('gene_metadata.db')
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS gene_info (
            gene_symbol TEXT PRIMARY KEY,
            name TEXT,
            type TEXT,
            aliases TEXT,
            summary TEXT,
            kegg_pathways TEXT,
            go_bp TEXT,
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    conn.commit()
    return conn

def fetch_batch(genes):
    url = "https://mygene.info/v3/query"
    payload = {
        'q': ','.join(genes),
        'scopes': 'symbol',
        'fields': 'name,symbol,type_of_gene,summary,alias,go.BP,pathway.kegg',
        'species': 'human',
        'dotfield': 'true'
    }
    
    try:
        r = requests.post(url, data=payload, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        print(f"Error fetching batch: {e}")
        return []

def main():
    print('Loading unique genes from midbase...')
    mid_conn = sqlite3.connect('midbase.db')
    genes_df = pd.read_sql('SELECT DISTINCT gene_symbol FROM differential_expression', mid_conn)
    all_genes = genes_df['gene_symbol'].tolist()
    
    meta_conn = setup_db()
    c = meta_conn.cursor()
    
    c.execute('SELECT gene_symbol FROM gene_info')
    existing = set(row[0] for row in c.fetchall())
    
    to_fetch = [g for g in all_genes if g not in existing]
    print(f"Total genes: {len(all_genes)}")
    print(f"Already cached: {len(existing)}")
    print(f"To fetch: {len(to_fetch)}")
    
    if not to_fetch:
        print("All genes cached!")
        return
        
    BATCH_SIZE = 1000
    
    print(f"Fetching in batches of {BATCH_SIZE}...")
    
    for i in range(0, len(to_fetch), BATCH_SIZE):
        batch = to_fetch[i:i+BATCH_SIZE]
        print(f"Processing batch {i//BATCH_SIZE + 1}/{(len(to_fetch)+BATCH_SIZE-1)//BATCH_SIZE} ({len(batch)} genes)")
        
        results = fetch_batch(batch)
        
        records_to_insert = []
        for hit in results:
            if 'notfound' in hit or 'symbol' not in hit:
                continue
                
            symbol = hit.get('symbol', '').upper()
            
            # Extract aliases
            aliases = hit.get('alias', [])
            if isinstance(aliases, str):
                aliases = [aliases]
                
            # Extract KEGG
            kegg_paths = []
            pathways = hit.get('pathway.kegg', [])
            if isinstance(pathways, dict):
                pathways = [pathways]
            for p in pathways:
                if isinstance(p, dict) and 'name' in p:
                    kegg_paths.append(p['name'])
                    
            # Extract GO BP
            go_terms = []
            go_bp = hit.get('go.BP', [])
            if isinstance(go_bp, dict):
                go_bp = [go_bp]
            for term in go_bp:
                if isinstance(term, dict) and 'term' in term:
                    go_terms.append(term['term'])
                    
            records_to_insert.append((
                symbol,
                hit.get('name', ''),
                hit.get('type_of_gene', ''),
                json.dumps(aliases),
                hit.get('summary', ''),
                json.dumps(kegg_paths),
                json.dumps(go_terms)
            ))
            
        if records_to_insert:
            c.executemany('''
                INSERT OR REPLACE INTO gene_info 
                (gene_symbol, name, type, aliases, summary, kegg_pathways, go_bp)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            ''', records_to_insert)
            meta_conn.commit()
            print(f"  Inserted {len(records_to_insert)} records")
            
        # Be nice to the API
        time.sleep(1)
        
    print("Done!")

if __name__ == "__main__":
    main()
