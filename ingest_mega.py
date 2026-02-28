import os
import pandas as pd
import sqlite3
import numpy as np
from scipy import stats
import json
import re

# Import the generated config
try:
    from src.data_config_full import FULL_NAMES_CONFIG
except ImportError:
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from src.data_config_full import FULL_NAMES_CONFIG

class MegaIngestionPipeline:
    def __init__(self, root_dir="data", db_path="midbase_core.db"):
        self.root_dir = root_dir
        self.processed_dir = f"{self.root_dir}/processed"
        self.matrix_path = f"{self.root_dir}/geo_downloads/matrix.tsv"
        self.metadata_path = f"{self.root_dir}/geo_downloads/Testis.json"
        self.db_path = db_path
        
        os.makedirs(self.processed_dir, exist_ok=True)
        
        # Load metadata once
        with open(self.metadata_path, 'r') as f:
            self.metadata = json.load(f)

    def load_master_matrix(self):
        print(f"📦 Loading Master Matrix from {self.matrix_path}...")
        try:
            df = pd.read_csv(self.matrix_path, sep='\t', index_col=0)
            print(f"   ✅ Loaded {df.shape[0]} genes x {df.shape[1]} samples")
            return df
        except Exception as e:
            print(f"   ❌ Error loading matrix: {e}")
            return None

    def parse_timepoint_name(self, gsm, gse_id):
        """Parse time-series metadata for meaningful naming."""
        if gse_id == "GSE116278" and gsm in self.metadata:
            # Fetal time series: parse age "6GW+0d" -> "Fetal_Week_6_Day_0"
            age_str = str(self.metadata[gsm].get('characteristics', ''))
            if 'age:' in age_str:
                age_val = age_str.split('age:')[-1].split(',')[0].strip()
                match = re.match(r'(\d+)GW\+?(\d*)d?', age_val)
                if match:
                    week = match.group(1)
                    day = match.group(2) if match.group(2) else '0'
                    return f"Fetal_Week_{week}_Day_{day}"
        return None

    def run(self):
        # Init DB
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute("DROP TABLE IF EXISTS studies")
        cursor.execute("CREATE TABLE studies (study_id TEXT PRIMARY KEY, title TEXT, category TEXT)")
        
        cursor.execute("DROP TABLE IF EXISTS samples")
        cursor.execute("CREATE TABLE samples (sample_id TEXT PRIMARY KEY, study_id TEXT, condition_group TEXT)")
        
        cursor.execute("DROP TABLE IF EXISTS differential_expression")
        cursor.execute("CREATE TABLE differential_expression (id INTEGER PRIMARY KEY, study_id TEXT, gene_symbol TEXT, comparison TEXT, logFC REAL, ave_expr REAL, p_value REAL)")
        
        conn.commit()

        # Load Data
        master_df = self.load_master_matrix()
        if master_df is None: return

        # Process Config
        for gse_id, info in FULL_NAMES_CONFIG.items():
            full_name = info['title']
            print(f"\n🚀 Processing {gse_id}: {full_name}...")
            
            # Identify Samples
            target_map = {}
            for group, gsms in info['samples'].items():
                for gsm in gsms: target_map[gsm] = group
            
            # Intersect
            valid_samples = [s for s in target_map.keys() if s in master_df.columns]
            
            if not valid_samples:
                print(f"   ⚠️ No samples found in matrix for {gse_id}. Skipping.")
                continue
                
            print(f"   Found {len(valid_samples)} samples.")
            
            # Extract & Transpose
            subset_df = master_df[valid_samples].T
            
            # Rename Samples
            new_index = []
            group_counts = {}
            original_gsms = subset_df.index.tolist()
            
            for gsm in original_gsms:
                group = target_map[gsm]
                
                # Try time-series naming first
                timepoint_name = self.parse_timepoint_name(gsm, gse_id)
                if timepoint_name:
                    new_id = timepoint_name
                else:
                    # Standard sequential naming
                    clean_group = group.replace(" ", "_").replace("/", "-")
                    if clean_group not in group_counts: group_counts[clean_group] = 0
                    group_counts[clean_group] += 1
                    new_id = f"{clean_group}_{group_counts[clean_group]}"
                
                new_index.append(new_id)
            
            subset_df.index = new_index
            subset_df['Condition_Label'] = [target_map[g] for g in original_gsms]
            
            # Export CSV
            csv_path = f"{self.processed_dir}/{gse_id}.csv"
            subset_df.to_csv(csv_path)
            print(f"   💾 Saved CSV: {csv_path} (Samples Renamed)")
            
            # DB Population
            cursor.execute("INSERT OR REPLACE INTO studies VALUES (?, ?, ?)", (gse_id, full_name, info['category']))
            sample_rows = [(gsm, gse_id, target_map[gsm]) for gsm in valid_samples]
            cursor.executemany("INSERT OR REPLACE INTO samples VALUES (?, ?, ?)", sample_rows)
            
            print(f"   ⏩ Skipped DE calculation (enable later if needed)")
            conn.commit()
            
        conn.close()
        print("\n✅ Mega-Ingestion Complete.")

if __name__ == "__main__":
    pipeline = MegaIngestionPipeline()
    pipeline.run()
