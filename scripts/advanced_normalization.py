"""
MIDBASE Advanced Co-Normalization Module
----------------------------------------
This script provides computationally intensive normalization techniques for 
projecting novel (N=1) patient RNA-seq samples into the MIDBASE reference space.

Due to the heavy memory requirements of these algorithms, they are excluded from 
the live Streamlit Community Cloud deployment and are provided here for local 
bioinformatics workflows and offline diagnostic pipelines.
"""

import pandas as pd
import numpy as np

def rank_based_normalization(reference_df, user_df):
    """
    Transforms both the reference database and the user sample into Percentile Ranks.
    
    Why use this? 
    Rank normalization completely ignores absolute read depth and platform-specific 
    background noise. A gene at the 99th percentile on a Microarray is directly 
    comparable to a gene at the 99th percentile on an Illumina NovaSeq.
    
    Parameters:
    -----------
    reference_df : pd.DataFrame
        The MIDBASE expression matrix (genes as index, samples as columns).
    user_df : pd.DataFrame
        The user's uploaded sample (genes as index, samples as columns).
        
    Returns:
    --------
    ranked_ref, ranked_user : tuple of pd.DataFrame
    """
    print("Aligning genes and applying Rank-Based Normalization...")
    
    # 1. Find common genes
    common_genes = reference_df.index.intersection(user_df.index)
    if len(common_genes) < 1000:
        raise ValueError("Insufficient common genes for reliable ranking.")
        
    ref_common = reference_df.loc[common_genes]
    user_common = user_df.loc[common_genes]
    
    # 2. Convert to percentile ranks (0.0 to 1.0) down each column
    ranked_ref = ref_common.rank(pct=True, method='average')
    ranked_user = user_common.rank(pct=True, method='average')
    
    print("Rank normalization complete.")
    return ranked_ref, ranked_user


def reference_quantile_normalization(reference_df, user_df):
    """
    Forces the user's uploaded sample to adopt the exact statistical distribution 
    (mean, median, variance curve) of the MIDBASE reference database.
    
    Why use this?
    This is the gold standard for projecting an N=1 sample into a historical PCA space.
    It removes technical batch effects while preserving the biological order of the user's genes.
    
    Parameters:
    -----------
    reference_df : pd.DataFrame
        The MIDBASE expression matrix.
    user_df : pd.DataFrame
        The user's uploaded sample.
        
    Returns:
    --------
    quantile_norm_user : pd.DataFrame
        The user's data, reshaped to fit the MIDBASE distribution.
    """
    print("Calculating MIDBASE reference distribution...")
    
    # 1. Find common genes
    common_genes = reference_df.index.intersection(user_df.index)
    ref_common = reference_df.loc[common_genes]
    user_common = user_df.loc[common_genes]
    
    # 2. Build the Target Distribution from the Reference
    # Sort every sample in the reference database lowest to highest
    sorted_ref = np.sort(ref_common.values, axis=0)
    # Take the mean across rows to get the "average" expression curve of the database
    target_distribution = np.mean(sorted_ref, axis=1)
    
    print("Applying Reference Quantile Normalization to user data...")
    # 3. Apply target distribution to the user's sample(s)
    q_norm_user = pd.DataFrame(index=user_common.index, columns=user_common.columns)
    
    for col in user_common.columns:
        sample_data = user_common[col]
        
        # Get the ranks of the user's genes (1 to N)
        ranks = sample_data.rank(method='min').astype(int) - 1
        
        # Map the target distribution values back to the user's genes based on their rank
        mapped_values = target_distribution[ranks]
        q_norm_user[col] = mapped_values
        
    print("Quantile normalization complete.")
    return q_norm_user


if __name__ == "__main__":
    # ==========================================
    # Example Usage for Local Deployment
    # ==========================================
    import os
    
    # Mock paths (Replace with actual paths when running)
    MIDBASE_MATRIX_PATH = "../data/processed/midbase_master_matrix.csv"
    USER_UPLOAD_PATH = "../sample_expression_data.csv"
    
    if os.path.exists(MIDBASE_MATRIX_PATH) and os.path.exists(USER_UPLOAD_PATH):
        # Load data (assuming genes are rows, samples are columns)
        ref_db = pd.read_csv(MIDBASE_MATRIX_PATH, index_col=0)
        user_data = pd.read_csv(USER_UPLOAD_PATH, index_col=0)
        
        # Method 1: Quantile Normalization (Best for PCA plotting)
        norm_user_q = reference_quantile_normalization(ref_db, user_data)
        
        # Method 2: Rank Normalization (Best for non-parametric machine learning)
        norm_ref_r, norm_user_r = rank_based_normalization(ref_db, user_data)
        
        # Save output for downstream PCA or ML models
        norm_user_q.to_csv("user_data_quantile_normalized.csv")
        print("Outputs saved successfully.")
    else:
        print("Run this script locally with valid data paths.")
