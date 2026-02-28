import requests
import json
import pandas as pd

def test_enrichr_api():
    gene_list = ['TP53', 'TNF', 'EGFR', 'AR']
    print(f"Testing Enrichr API with {len(gene_list)} genes: {gene_list}")
    
    try:
        # 1. Add list to Enrichr
        add_url = "https://maayanlab.cloud/Enrichr/addList"
        # The 'list' part of the payload needs to be a string of newline-separated genes.
        # It's passed as a 'file' parameter in requests.
        payload = {
            'list': (None, '\n'.join(gene_list)),
            'description': (None, 'KG_Network_Genes')
        }
        print(f"Sending POST to {add_url}")
        res = requests.post(add_url, files=payload)
        print(f"POST status: {res.status_code}")
        
        if not res.ok: 
            print(f"POST failed: {res.text}")
            return
            
        data = res.json()
        print(f"POST response JSON: {data}")
        user_list_id = data['userListId']
        
        # 2. Get enrichment results
        libraries = ['KEGG_2021_Human', 'GO_Biological_Process_2023']
        
        for lib in libraries:
            enrich_url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={lib}"
            print(f"\nSending GET to {enrich_url}")
            res_enrich = requests.get(enrich_url)
            print(f"GET status: {res_enrich.status_code}")
            
            if not res_enrich.ok: 
                print(f"GET failed: {res_enrich.text}")
                continue
                
            enrich_data = res_enrich.json()
            if lib not in enrich_data:
                print(f"Library '{lib}' not in response JSON keys: {list(enrich_data.keys())}")
                continue
                
            results = enrich_data[lib]
            print(f"Found {len(results)} raw results for {lib}.")
            if results:
                print(f"First result sample: {results[0]}")
            
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    test_enrichr_api()
