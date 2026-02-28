import json
import os

# User-Defined Full Descriptive Names & Specs
DATASETS = {
    # Disease / Conditions
    "GSE106487": {
        "name": "Obstructive_Azoospermia_vs_Non_Obstructive_Azoospermia",
        "groups": {
            "Obstructive_Azoospermia": 61,
            "Non_Obstructive_Azoospermia": 4
        },
        "category": "Disease"
    },
    "GSE103905": {
        "name": "Klinefelter_Syndrome_and_Sertoli_Cell_Only_Syndrome",
        "groups": {
            "Normal_Spermatogenesis": 8,
            "Klinefelter_Syndrome": 10,
            "Sertoli_Cell_Only_Syndrome": 4
        },
        "category": "Disease"
    },
    "GSE154535": {
        "name": "Idiopathic_Non_Obstructive_Azoospermia_vs_Obstructive",
        "groups": {
            "Obstructive_Azoospermia": 1,
            "Idiopathic_Non_Obstructive_Azoospermia": 2
        },
        "category": "Disease"
    },
    "GSE125222": {
        "name": "Androgen_Insensitivity_Syndrome",
        "groups": {
            "Normal_Control": 3,
            "Androgen_Insensitivity_Syndrome": 3
        },
        "category": "Disease"
    },
    "GSE200680": {
        "name": "Testicular_Fibrosis",
        "groups": {
            "Non_Fibrotic_Control": 10,
            "Fibrotic_Testicular_Tissue": 10
        },
        "category": "Disease"
    },
    "GSE208761": {
        "name": "Klinefelter_Syndrome_Cultured_Cells",
        "groups": {
            "XY_Control": 1,
            "XXY_Klinefelter": 1
        },
        "category": "Disease"
    },
    "GSE235210": {
        "name": "Zika_Virus_Infection_Model",
        "groups": {
            "Mock_Treated_Control": 3,
            "Zika_Virus_Infection": 12
        },
        "category": "Disease"
    },
    
    # Development
    "GSE116278": {
        "name": "Fetal_Testis_Development_Time_Series",
        "groups": {"Fetal_Timeseries": "all"}, # All fetal samples (6-17 weeks)
        "category": "Development"
    },
    "GSE124263": {
        "name": "Neonatal_vs_Adult_Testis",
        "groups": {
            "Neonatal": "heuristic",
            "Adult": "heuristic"
        },
        "category": "Development"
    },
    
    # Stem Cells / Cell Types
    "GSE108977": {
        "name": "Spermatogonial_Stem_Cells_ID4_GFP",
        "groups": {"Spermatogonia": 635},
        "category": "Cell Type"
    },
    "GSE92276": {
        "name": "Spermatogonial_Stem_Cells_Sorted",
        "groups": {"Spermatogonial_Stem_Cells": 368},
        "category": "Cell Type"
    },
    # ... Add others as we verify they exist
}

def generate_config():
    base_dir = "d:/Projects/Merge_Midbase_Serenova/XY_Counsel/data/geo_downloads"
    json_path = os.path.join(base_dir, "Testis.json")
    
    print(f"📂 Loading {json_path}...")
    with open(json_path, 'r') as f:
        data = json.load(f)
        
    print(f"   Total Samples in Metadata: {len(data)}")
    
    generated_config = {}
    
    # Index by Series
    series_map = {} # GSE -> [list of sample dicts]
    for gsm, info in data.items():
        # Info structure check: "series" is string or list?
        if 'series' not in info: continue
        series = info['series']
        # Handle list vs string
        series_list = series if isinstance(series, list) else [series]
        # remove comma split strings just in case "GSE1,GSE2"
        all_series = []
        for s in series_list:
            all_series.extend(s.split(','))
            
        for s in all_series:
            s = s.strip()
            if s not in series_map: series_map[s] = []
            info['gsm'] = gsm
            series_map[s].append(info)

    # Process each requested dataset
    for gse, spec in DATASETS.items():
        print(f"\n🔍 Processing {gse} ({spec['name']})...")
        if gse not in series_map:
            print(f"   ⚠️ {gse} NOT FOUND in metadata!")
            continue
            
        samples = series_map[gse]
        print(f"   Found {len(samples)} samples metadata.")
        
        # Keyword matching / Grouping logic
        mapped_samples = {}
        
        # Strategy: Iterate defined groups and try to find them
        # If group count is integer, try to find that many matching checks
        # If group is "all", take all
        
        used_gsms = set()
        
        for group_name, count_spec in spec['groups'].items():
            found_gsms = []
            
            # Simple keyword match on 'characteristics' or 'source'
            keywords = group_name.replace('_', ' ').split()
            
            # Filter samples
            candidates = []
            for s in samples:
                if s['gsm'] in used_gsms: continue
                
                # Create searchable text
                text = str(s.get('characteristics', '')) + " " + str(s.get('source_name_ch1', '')) + " " + str(s.get('title', ''))
                text = text.lower()
                
                # Heuristic: Check if keywords match
                # For "Normal Control", look for "normal" or "control"
                # For "Klinefelter", look for "klinefelter" or "xxy"
                
                match_score = 0
                g_lower = group_name.lower()
                
                if "klinefelter" in g_lower or "xxy" in g_lower:
                    if "klinefelter" in text or "xxy" in text or "47,xxy" in text: match_score += 1
                elif "sertoli" in g_lower:
                    if "sertoli" in text and "only" in text: match_score += 1
                elif "fibrosis" in g_lower:
                    if "fibro" in text: match_score += 1
                elif "spermatogonial" in g_lower:
                    # GSE92276, GSE108977: Look for "SSC", "Spermatogonia", "Spermatogonium", "SSEA4", "ID4", "KIT"
                    if "spermatogonia" in text or "spermatogonium" in text or "ssc" in text or "ssea" in text or "id4" in text or "kit" in text: match_score += 1
                elif "obstructive" in g_lower and "non" not in g_lower:
                    # OA (Exclude NOA) - Look for "OA" or "obstructive" WITHOUT "non"
                    if ("oa" in text and "noa" not in text) or ("obstructive" in text and "non" not in text):
                        match_score += 1
                elif "idiopathic_non_obstructive" in g_lower:
                     # GSE154535 Specific
                     if "idiopathic" in text and ("noa" in text or "non-obstructive" in text):
                         match_score += 1
                elif "non_obstructive" in g_lower:
                    # NOA Generic - Look for "NOA" or "non-obstructive" or "non obstructive"
                    if "noa" in text or "non-obstructive" in text or "non obstructive" in text:
                        match_score += 1
                elif "androgen" in g_lower and "insensitivity" in g_lower:
                    # AIS - Look for "androgen" + "insensitivity" or "ais"
                    if ("androgen" in text and "insensitivity" in text) or "ais" in text:
                        match_score += 1
                elif "normal" in g_lower or "control" in g_lower or "xy_control" in g_lower:
                    if "normal" in text or "control" in text or "healthy" in text or ("xy" in text and "xxy" not in text):
                        match_score += 1
                elif "zika" in g_lower:
                    if "zika" in text or "zikv" in text: match_score += 1
                elif "mock" in g_lower:
                    if "mock" in text: match_score += 1
                elif "fetal" in g_lower or "timeseries" in g_lower:
                    # GSE116278 - Look for "GW" (gestational weeks), "fetal", "embryo"
                    if "gw" in text or "gestational" in text or "fetal" in text or "embryo" in text:
                        match_score += 1
                elif "prenatal" in g_lower:
                    if "prenatal" in text or "week" in text: match_score += 1
                elif "neonatal" in g_lower:
                    if "neonatal" in text or "infant" in text or "day" in text: match_score += 1
                elif "adult" in g_lower:
                    if "adult" in text or "year" in text: match_score += 1
                
                if match_score > 0:
                    candidates.append(s['gsm'])
            
            # DEBUG: If 0 found, print first few samples metadata to help debug
            if len(candidates) == 0:
                 print(f"      ❌ Zero matches for '{group_name}'. Inspecting first 3 samples of {gse}:")
                 for i, s in enumerate(samples[:3]):
                      t = str(s.get('title', ''))
                      src = str(s.get('source_name_ch1', ''))
                      char = str(s.get('characteristics', ''))
                      print(f"         [{i}] Title: {t} | Source: {src} | Char: {char}")

            # If strictly N samples requested, take top N or all if candidates
            if isinstance(count_spec, int):
                # Try to fit exact
                if len(candidates) >= count_spec:
                    found_gsms = candidates[:count_spec]
                else:
                    print(f"      Running wide search for {group_name} (Need {count_spec}, found {len(candidates)})")
                    # Fallback logic could go here
                    found_gsms = candidates
            elif count_spec == "all":
                found_gsms = candidates
                
            print(f"   -> {group_name}: Found {len(found_gsms)} (Target: {count_spec})")
            mapped_samples[group_name] = found_gsms
            for g in found_gsms: used_gsms.add(g)
            
        generated_config[gse] = {
            "title": spec['name'], # Full Name
            "samples": mapped_samples,
            "category": spec['category']
        }

    # Output to python file
    out_path = "src/data_config_full.py"
    with open(out_path, 'w') as f:
        f.write("# Generated Data Config with Full Descriptive Names\n\n")
        f.write("FULL_NAMES_CONFIG = ")
        f.write(json.dumps(generated_config, indent=4))
        
    print(f"\n✅ Generated config at {out_path}")

if __name__ == "__main__":
    generate_config()
