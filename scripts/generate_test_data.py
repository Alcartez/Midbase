import pandas as pd
import numpy as np

# Load real midbase data to use as a baseline
try:
    # The instruction implies a db_path variable should exist and be set to "midbase_core.db"
    # Although not explicitly used in the subsequent CSV loading, this line is added as per instruction.
    db_path = "midbase_core.db"

    klinefelter_data = pd.read_csv("data/processed/GSE103905.csv", index_col=0)
    noa_data = pd.read_csv("data/processed/GSE106487.csv", index_col=0)
    
    # Extract one sample of each condition type, ensuring we select actual disease samples
    # In these processed CSVs, samples are rows (index) and genes are columns
    klinefelter_sample = klinefelter_data.loc[["Klinefelter_Syndrome_1"]].T
    noa_sample = noa_data.loc[["Non_Obstructive_Azoospermia_1"]].T
    
    # Merge into a single "User Uploaded" dataframe
    # We use an inner join to ensure the genes exist in both
    user_test_df = pd.merge(klinefelter_sample, noa_sample, left_index=True, right_index=True)
    
    # Drop condition labels if they were accidentally transposed in as a row
    if 'Condition_Label' in user_test_df.index:
        user_test_df = user_test_df.drop('Condition_Label')
        
    user_test_df = user_test_df.astype(float)
    
    # Rename columns to simulate an external lab uploading anonymous patients
    user_test_df.columns = ["Patient_A_Suspected_Klinefelter", "Patient_B_Suspected_NOA"]
    
    # Add varying degrees of random multiplicative noise to simulate sequencing depth differences
    # Multiplicative noise ensures genes that are strictly 0 (not expressed) stay 0,
    # preventing artificial clustering in PCA due to thousands of unexpressed genes gaining noise variance.
    noise_A = np.random.normal(1.0, 0.05, len(user_test_df))
    noise_B = np.random.normal(1.0, 0.08, len(user_test_df))
    
    user_test_df["Patient_A_Suspected_Klinefelter"] *= noise_A
    user_test_df["Patient_B_Suspected_NOA"] *= noise_B
    
    # Ensure no negative expression values after noise
    user_test_df[user_test_df < 0] = 0
    
    # Export to CSV for the user to upload
    output_path = "sample_expression_data.csv"
    user_test_df.to_csv(output_path)
    print(f"✅ Successfully generated realistic test data: {output_path}")
    print(f"Shape: {user_test_df.shape[0]} genes, {user_test_df.shape[1]} samples")

except Exception as e:
    print(f"Error generating test data: {e}")
