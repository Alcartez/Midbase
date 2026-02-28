"""
Create Simple Test File
"""
import pandas as pd
import numpy as np

# Load a sample from GSE103905 (genes should already be columns)
df = pd.read_csv("data/processed/GSE103905.csv", index_col=0)

# Drop Condition_Label
if 'Condition_Label' in df.columns:
    df = df.drop(columns=['Condition_Label'])

print(f"Loaded data shape: {df.shape}")
print(f"Rows (samples): {len(df)}, Columns (genes): {len(df.columns)}")

# Select one sample (row) and create noisy copies
template = df.iloc[0]  # First sample

# Create 2 noisy versions
np.random.seed(42)
sample1 = template + np.random.normal(0, template.std() * 0.05, len(template))
sample2 = template + np.random.normal(0, template.std() * 0.08, len(template))

# Ensure non-negative  
sample1 = np.maximum(sample1, 0.01)
sample2 = np.maximum(sample2, 0.01)

# Create new dataframe (samples as rows, genes as columns)
test_df = pd.DataFrame([sample1, sample2], columns=df.columns)
test_df.index = ['TestSample_1', 'TestSample_2']

# Transpose so genes are rows, samples are columns (expected upload format)
test_df = test_df.T

# Save
test_df.to_csv("sample_expression_data.csv")

print(f"\n✅ Created sample_expression_data.csv")
print(f"   Shape: {test_df.shape[0]} genes × {test_df.shape[1]} samples")
print(f"   Mean expression: {test_df.mean().mean():.2f}")
print(f"   First few genes: {list(test_df.index[:5])}")
print("\n📤 Upload this file to test the Upload & Compare feature!")
print("   Expected match: GSE103905 (Klinefelter study)")

