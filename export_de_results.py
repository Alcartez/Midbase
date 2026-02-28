"""
Export Differential Expression Results to CSV
Exports DE statistics from midbase.db to individual CSV files per study
and a combined master CSV with all results.
"""
import os
import sqlite3
import pandas as pd

class DEExporter:
    def __init__(self, db_path="midbase_de.db",
                 output_dir="de_results"):
        self.db_path = db_path
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def export_all(self):
        """Export all DE results from database."""
        print(f"📊 Connecting to database: {self.db_path}")
        conn = sqlite3.connect(self.db_path)
        
        # Load all DE results
        query = """
        SELECT 
            de.study_id,
            s.title as study_title,
            de.gene_symbol,
            de.logFC,
            de.ave_expr,
            de.p_value,
            de.comparison
        FROM differential_expression de
        LEFT JOIN studies s ON de.study_id = s.study_id
        ORDER BY de.study_id, de.p_value
        """
        
        print("📥 Loading differential expression results...")
        df_all = pd.read_sql_query(query, conn)
        print(f"   ✅ Loaded {len(df_all)} DE results across {df_all['study_id'].nunique()} studies")
        
        # Export combined results
        combined_path = os.path.join(self.output_dir, "All_DE_Results.csv")
        df_all.to_csv(combined_path, index=False)
        print(f"\n💾 Saved combined results: {combined_path}")
        
        # Export individual study files
        print(f"\n📁 Exporting individual study files...")
        for study_id in df_all['study_id'].unique():
            study_df = df_all[df_all['study_id'] == study_id].copy()
            
            # Add adjusted p-value (Benjamini-Hochberg FDR)
            study_df = study_df.sort_values('p_value')
            n = len(study_df)
            study_df['rank'] = range(1, n + 1)
            study_df['adj_p_value'] = study_df['p_value'] * n / study_df['rank']
            study_df['adj_p_value'] = study_df['adj_p_value'].clip(upper=1.0)
            
            # Sort by significance
            study_df = study_df.sort_values('adj_p_value')
            study_df = study_df.drop(columns=['rank'])
            
            # Save
            study_path = os.path.join(self.output_dir, f"{study_id}_DE_Results.csv")
            study_df.to_csv(study_path, index=False)
            
            sig_count = (study_df['adj_p_value'] < 0.05).sum()
            print(f"   {study_id}: {len(study_df)} genes ({sig_count} significant at FDR < 0.05)")
        
        # Export top genes summary
        print(f"\n🔝 Creating top differentially expressed genes summary...")
        top_genes = []
        for study_id in df_all['study_id'].unique():
            study_df = df_all[df_all['study_id'] == study_id].copy()
            
            # Add adjusted p-values
            study_df = study_df.sort_values('p_value')
            n = len(study_df)
            study_df['rank'] = range(1, n + 1)
            study_df['adj_p_value'] = study_df['p_value'] * n / study_df['rank']
            study_df['adj_p_value'] = study_df['adj_p_value'].clip(upper=1.0)
            
            # Filter significant and top by absolute logFC
            sig_df = study_df[study_df['adj_p_value'] < 0.05].copy()
            if len(sig_df) > 0:
                sig_df['abs_logFC'] = sig_df['logFC'].abs()
                top_20 = sig_df.nlargest(20, 'abs_logFC')
                top_genes.append(top_20[['study_id', 'study_title', 'gene_symbol', 'logFC', 'p_value', 'adj_p_value']])
        
        if top_genes:
            top_df = pd.concat(top_genes, ignore_index=True)
            top_path = os.path.join(self.output_dir, "Top_DE_Genes_Summary.csv")
            top_df.to_csv(top_path, index=False)
            print(f"   💾 Saved top genes summary: {top_path}")
            print(f"   📊 Total significant genes across all studies: {len(top_df)}")
        
        conn.close()
        print(f"\n✅ Export complete! All files saved to: {self.output_dir}")
        
        # Print summary statistics
        print(f"\n📈 Summary Statistics:")
        print(f"   Total DE comparisons: {len(df_all)}")
        print(f"   Studies with DE results: {df_all['study_id'].nunique()}")
        print(f"   Unique genes analyzed: {df_all['gene_symbol'].nunique()}")

if __name__ == "__main__":
    exporter = DEExporter()
    exporter.export_all()
