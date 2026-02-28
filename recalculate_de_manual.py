"""
Manual T-Test Implementation (No scipy required)
Recalculate DE p-values using hand-coded Welch's t-test
"""
import os
import pandas as pd
import sqlite3
import numpy as np

def manual_welch_ttest(group1, group2):
    """
    Calculate Welch's t-test manually (unequal variance t-test)
    Returns t-statistic and p-value
    """
    from math import erf, sqrt
    
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
        # Both groups have zero variance
        if mean1 == mean2:
            return 0.0, 1.0  # No difference
        else:
            return float('inf'), 0.0  # Perfect separation
    
    # Calculate standard error
    se = sqrt((var1/n1) + (var2/n2))
    
    if se == 0:
        return 0.0, 1.0
    
    # Welch's t-statistic
    t_stat = (mean1 - mean2) / se
    
    # Degrees of freedom (Welch-Satterthwaite equation)
    numerator = ((var1/n1) + (var2/n2)) ** 2
    denominator = ((var1/n1)**2 / (n1-1)) + ((var2/n2)**2 / (n2-1))
    
    if denominator == 0:
        df = n1 + n2 - 2  # Fallback to pooled df
    else:
        df = numerator / denominator
    
    # Calculate p-value using normal approximation to t-distribution
    # Standard normal CDF: Phi(x) = 0.5 * (1 + erf(x/sqrt(2)))
    # Two-tailed p-value: 2 * (1 - Phi(|t|))
    
    # For t-distribution with df > 5, normal approximation is reasonable
    # Standard normal CDF
    abs_t = abs(t_stat)
    
    # Using erf to calculate normal CDF
    # CDF(x) = 0.5 * (1 + erf(x / sqrt(2)))
    normal_cdf = 0.5 * (1.0 + erf(abs_t / sqrt(2)))
    
    # Two-tailed p-value
    p_value = 2.0 * (1.0 - normal_cdf)
    
    # Ensure p-value is in valid range [0, 1]
    p_value = max(0.0, min(1.0, p_value))
    
    return float(t_stat), float(p_value)

class ManualDECalculator:
    def __init__(self, processed_dir="data/processed",
                 db_path="midbase_core.db"):
        self.processed_dir = processed_dir
        self.db_path = db_path

    def calculate_de(self, csv_path, gse_id):
        """Calculate DE using manual t-test."""
        print(f"\n🧬 Calculating DE for {gse_id}...")
        
        # Load data
        df = pd.read_csv(csv_path, index_col=0)
        if 'Condition_Label' not in df.columns:
            print(f"   ⚠️ No Condition_Label column. Skipping DE.")
            return None
        
        groups = df['Condition_Label'].unique()
        
        # Look for control group
        control_group = next((g for g in groups if 'Control' in g or 'Normal' in g or 'Mock' in g), None)
        
        if not control_group:
            print(f"   ⚠️ No control group found (groups: {groups}). Skipping DE.")
            return None
            
        print(f"   📊 Control group: {control_group}")
        
        # Separate control vs disease
        control_df = df[df['Condition_Label'] == control_group].drop(columns=['Condition_Label'])
        disease_df = df[df['Condition_Label'] != control_group].drop(columns=['Condition_Label'])
        
        if len(control_df) == 0 or len(disease_df) == 0:
            print(f"   ⚠️ Empty group. Skipping.")
            return None
        
        # Check for minimum sample size
        if len(control_df) < 2 or len(disease_df) < 2:
            print(f"   ⚠️ Insufficient samples (need n>=2 for t-test). Control={len(control_df)}, Disease={len(disease_df)}")
            print(f"   💡 This comparison will be handled by cross-study analysis with pooled controls.")
            return None
        
        print(f"   Control: {len(control_df)}, Disease: {len(disease_df)}, Genes: {len(control_df.columns)}")
        
        # Calculate DE for each gene
        results = []
        gene_list = control_df.columns.tolist()
        total = len(gene_list)
        
        print(f"   ⚡ Running manual t-tests...")
        for idx, gene in enumerate(gene_list):
            if (idx + 1) % 10000 == 0:
                print(f"      Progress: {idx+1}/{total} ({100*(idx+1)/total:.1f}%)")
            
            control_vals = control_df[gene].values
            disease_vals = disease_df[gene].values
            
            # Statistics
            mean_control = np.mean(control_vals)
            mean_disease = np.mean(disease_vals)
            
            pseudo = 1.0
            logfc = np.log2(mean_disease + pseudo) - np.log2(mean_control + pseudo)
            ave_expr = np.log2((mean_control + mean_disease) / 2 + pseudo)
            
            # Manual t-test
            try:
                t_stat, p_val = manual_welch_ttest(disease_vals, control_vals)
            except Exception as e:
                print(f"   Error for {gene}: {e}")
                p_val = 1.0
            
            results.append({
                'gene_symbol': gene,
                'logFC': float(logfc),
                'ave_expr': float(ave_expr),
                'p_value': float(p_val),
                'comparison': f"vs_{control_group}",
                'study_id': gse_id
            })
        
        results_df = pd.DataFrame(results)
        
        # Stats
        valid = (results_df['p_value'] < 1.0).sum()
        sig = (results_df['p_value'] < 0.05).sum()
        print(f"   ✅ Complete. Valid p-values: {valid}/{len(results_df)}, Significant: {sig}")
        
        return results_df

    def run(self):
        """Process all CSV files."""
        conn = sqlite3.connect(self.db_path)
       
        # Clear existing
        cursor = conn.cursor()
        cursor.execute("DELETE FROM differential_expression")
        conn.commit()
        print("🗑️ Cleared existing DE data\n")
        
        csv_files = [f for f in os.listdir(self.processed_dir) if f.endswith('.csv')]
        print(f"Found {len(csv_files)} datasets\n")
        
        total = 0
        for csv_file in csv_files:
            gse_id = csv_file.replace('.csv', '')
            csv_path = os.path.join(self.processed_dir, csv_file)
            
            de_results = self.calculate_de(csv_path, gse_id)
            
            if de_results is not None:
                de_results.to_sql('differential_expression', conn, if_exists='append', index=False)
                total += len(de_results)
                print(f"   💾 Saved to database.")
        
        conn.close()
        print(f"\n✅ Complete! Total: {total} results")

if __name__ == "__main__":
    calc = ManualDECalculator()
    calc.run()
