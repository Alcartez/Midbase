"""
Generate Test Data for Upload & Compare Feature
Creates sample gene expression data for testing the upload functionality
"""
import pandas as pd
import numpy as np
import sqlite3

def generate_test_data(output_file="test_expression_data.csv", n_samples=3, similarity_study='GSE103905'):
    """
    Generate test expression data that resembles a specific study
    
    Parameters:
    - output_file: CSV filename to save
    - n_samples: Number of samples to generate
    - similarity_study: Which study to base the expression on
    """
    
    # Load midbase expression for the target study
    db_path = "d:/Projects/Merge_Midbase_Serenova/XY_Counsel/midbase.db"
    processed_dir = "d:/Projects/Merge_Midbase_Serenova/XY_Counsel/data/processed"
    
    print(f"Generating {n_samples} test samples similar to {similarity_study}...")
    
    # Load target study
    study_csv = f"{processed_dir}/{similarity_study}.csv"
    df = pd.read_csv(study_csv, index_col=0)
    
    # Drop Condition_Label if exists
    if 'Condition_Label' in df.columns:
        df = df.drop(columns=['Condition_Label'])
    
    # Transpose (so genes are columns, samples are rows)
    df = df.T
    
    # Select a random sample from the study as template
    template_sample = df.sample(1)
    
    # Generate new samples with noise
    test_samples = []
    for i in range(n_samples):
        # Add Gaussian noise to template
        noise_scale = 0.1  # 10% noise
        noisy_sample = template_sample.values[0] + np.random.normal(0, template_sample.values[0].std() * noise_scale, len(template_sample.columns))
        
        # Ensure non-negative
        noisy_sample = np.maximum(noisy_sample, 0.1)
        
        test_samples.append(noisy_sample)
    
    # Create DataFrame
    test_df = pd.DataFrame(test_samples, columns=template_sample.columns)
    test_df.index = [f'TestSample_{i+1}' for i in range(n_samples)]
    
    # Transpose back (genes as rows, samples as columns)
    test_df = test_df.T
    
    # Save to CSV
    test_df.to_csv(output_file)
    
    print(f"✅ Saved {len(test_df)} genes × {len(test_df.columns)} samples to {output_file}")
    print(f"   Data shape: {test_df.shape}")
    print(f"   Mean expression: {test_df.mean().mean():.2f}")
    print(f"   Based on study: {similarity_study}")
    print(f"\n📊 This data should match highly with {similarity_study} samples!")
    
    return test_df

if __name__ == "__main__":
    # Generate test data similar to Klinefelter study
    generate_test_data("test_klinefelter_like.csv", n_samples=2, similarity_study='GSE103905')
    
    # Generate test data similar to AIS study
    generate_test_data("test_ais_like.csv", n_samples=1, similarity_study='GSE125222')
    
    print("\n✅ Test files generated! Upload these to the platform to test the comparison feature.")
