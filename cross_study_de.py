"""
Cross-Study Differential Expression Analysis
Uses pooled normal controls from multiple studies to compare against disease samples
from studies lacking their own controls.
"""
import os
import pandas as pd
import sqlite3
import numpy as np
from math import erf, sqrt

def manual_welch_ttest(group1, group2):
    """
    Calculate Welch's t-test manually (unequal variance t-test)
    Returns t-statistic and p-value
    """
    n1 = len(group1)
    n2 = len(group2)
    
    # Sample means
    mean1 = np.mean(group1)
    mean2 = np.mean(group2)
    
    # Sample variances (unbiased estimator)
    var1 = np.var(group1, ddof=1)
    var2 = np.var(group2, ddof=1)
    
    # Handle edge cases
    if var1 == 0 and var2 == 0:
        if mean1 == mean2:
            return 0.0, 1.0
        else:
            return float('inf'), 0.0
    
    # Calculate standard error
    se = sqrt((var1/n1) + (var2/n2))
    
    if se == 0:
        return 0.0, 1.0
    
    # Welch's t-statistic
    t_stat = (mean1 - mean2) / se
    
    # Calculate p-value using normal approximation
    abs_t = abs(t_stat)
    normal_cdf = 0.5 * (1.0 + erf(abs_t / sqrt(2)))
    p_value = 2.0 * (1.0 - normal_cdf)
    p_value = max(0.0, min(1.0, p_value))
    
    return float(t_stat), float(p_value)

class CrossStudyDE:
    def __init__(self, processed_dir="data/processed",
                 db_path="midbase_core.db"):
        self.processed_dir = processed_dir
        self.db_path = db_path
        
        # Define control sources and disease targets
        self.control_sources = {
            'GSE103905': 'Normal_Spermatogenesis',
            'GSE125222': 'Normal_Control',
            'GSE208761': 'XY_Control'
        }
        
        self.disease_targets = {
            'GSE106487': {
                'groups': ['Obstructive_Azoospermia', 'Non_Obstructive_Azoospermia'],
                'comparisons': [
                    ('Pooled_Normal', 'Obstructive_Azoospermia'),
                    ('Pooled_Normal', 'Non_Obstructive_Azoospermia')
                ]
            },
            'GSE154535': {
                'groups': ['Obstructive_Azoospermia', 'Idiopathic_Non_Obstructive_Azoospermia'],
                'comparisons': [
                    ('Pooled_Normal', 'Obstructive_Azoospermia'),
                    ('Pooled_Normal', 'Idiopathic_Non_Obstructive_Azoospermia')
                ]
            },
            'GSE235210': {
                'groups': ['Zika_Virus_Infection', 'Mock_Treated_Control'],
                'comparisons': [
                    ('Pooled_Normal', 'Zika_Virus_Infection')
                ]
            }
        }

    def load_and_pool_controls(self):
        """Load and pool normal control samples from multiple studies."""
        print("\n📦 Loading normal control samples from multiple studies...")
        
        all_controls = []
        control_info = []
        
        for study_id, control_label in self.control_sources.items():
            csv_path = os.path.join(self.processed_dir, f"{study_id}.csv")
            if not os.path.exists(csv_path):
                print(f"   ⚠️ {study_id} not found. Skipping.")
                continue
            
            df = pd.read_csv(csv_path, index_col=0)
            if 'Condition_Label' not in df.columns:
                continue
            
            # Extract controls
            control_df = df[df['Condition_Label'] == control_label].drop(columns=['Condition_Label'])
            
            if len(control_df) > 0:
                all_controls.append(control_df)
                control_info.append({
                    'study': study_id,
                    'label': control_label,
                    'count': len(control_df)
                })
                print(f"   ✅ {study_id}: {len(control_df)} {control_label} samples")
        
        if not all_controls:
            print("   ❌ No control samples found!")
            return None, None
        
        # Combine all controls
        pooled_controls = pd.concat(all_controls, axis=0)
        print(f"\n   📊 Pooled Controls: {len(pooled_controls)} samples from {len(control_info)} studies")
        
        return pooled_controls, control_info

    def calculate_cross_study_de(self, control_df, disease_df, comparison_name, study_id):
        """Calculate DE between pooled controls and disease samples using manual t-test."""
        print(f"\n   🧬 {comparison_name}")
        print(f"      Control: {len(control_df)} samples, Disease: {len(disease_df)} samples")
        
        # Calculate DE for each gene
        results = []
        gene_list = control_df.columns.tolist()
        total = len(gene_list)
        
        print(f"      ⚡ Running manual t-tests...")
        for idx, gene in enumerate(gene_list):
            if (idx + 1) % 10000 == 0:
                print(f"         Progress: {idx+1}/{total} ({100*(idx+1)/total:.1f}%)")
            
            control_vals = control_df[gene].values
            disease_vals = disease_df[gene].values
            
            # Calculate statistics
            mean_control = np.mean(control_vals)
            mean_disease = np.mean(disease_vals)
            
            # Log fold change
            pseudo = 1.0
            logfc = np.log2(mean_disease + pseudo) - np.log2(mean_control + pseudo)
            ave_expr = np.log2((mean_control + mean_disease) / 2 + pseudo)
            
            # Manual t-test
            try:
                t_stat, p_val = manual_welch_ttest(disease_vals, control_vals)
            except Exception as e:
                p_val = 1.0
            
            results.append({
                'gene_symbol': gene,
                'logFC': float(logfc),
                'ave_expr': float(ave_expr),
                'p_value': p_val,
                'comparison': comparison_name,
                'study_id': study_id
            })
        
        results_df = pd.DataFrame(results)
        
        # Stats
        valid = (results_df['p_value'] < 1.0).sum()
        sig = (results_df['p_value'] < 0.05).sum()
        print(f"      ✅ Complete. Valid p-values: {valid}/{len(results_df)}, Significant: {sig}")
        
        return results_df

    def run(self):
        """Perform cross-study DE analysis."""
        print("🔬 Cross-Study Differential Expression Analysis")
        print("=" * 60)
        
        # Load pooled controls
        pooled_controls, control_info = self.load_and_pool_controls()
        if pooled_controls is None:
            return
        
        # Connect to database
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        all_results = []
        
        # Process each disease study
        for study_id, config in self.disease_targets.items():
            print(f"\n\n{'='*60}")
            print(f"📊 Processing {study_id}")
            print(f"{'='*60}")
            
            csv_path = os.path.join(self.processed_dir, f"{study_id}.csv")
            if not os.path.exists(csv_path):
                print(f"   ⚠️ File not found: {csv_path}")
                continue
            
            df = pd.read_csv(csv_path, index_col=0)
            if 'Condition_Label' not in df.columns:
                print(f"   ⚠️ No Condition_Label column")
                continue
            
            # For each comparison
            for control_label, disease_label in config['comparisons']:
                disease_df = df[df['Condition_Label'] == disease_label].drop(columns=['Condition_Label'])
                
                if len(disease_df) == 0:
                    print(f"   ⚠️ No samples for {disease_label}")
                    continue
                
                # Calculate DE
                comparison_name = f"Pooled_Normal_vs_{disease_label}"
                de_results = self.calculate_cross_study_de(
                    pooled_controls, 
                    disease_df, 
                    comparison_name, 
                    study_id
                )
                
                all_results.append(de_results)
                
                # Save to database
                de_results.to_sql('differential_expression', conn, if_exists='append', index=False)
                print(f"      💾 Saved to database")
        
        conn.close()
        
        if all_results:
            print(f"\n\n{'='*60}")
            print(f"✅ Cross-Study DE Analysis Complete!")
            print(f"{'='*60}")
            print(f"   Total comparisons: {len(all_results)}")
            print(f"   Total results: {sum(len(r) for r in all_results)}")
        else:
            print("\n⚠️ No results generated")

if __name__ == "__main__":
    analyzer = CrossStudyDE()
    analyzer.run()
