import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import glob


# Set input and output folders
input_folder = "Jaccardinput"  # Path to folder containing TSV files
output_folder = "jaccard_heatmaps_triangle"
binary_output_folder = "binary_matrices" #to get txt files with matrix of detected or not
os.makedirs(output_folder, exist_ok=True)
os.makedirs(binary_output_folder, exist_ok=True)

# Get all .tsv files in the folder
tsv_files = glob.glob(os.path.join(input_folder, "*_long.tsv"))

def compute_jaccard(df, file_name):
    # Create a unique isoform_id (Gene|Isoform|SampleType) to handle duplicate isoforms across sample types
    df['isoform_id'] = df['Gene'].astype(str) + "|" + df['Isoform'].astype(str) + "|" + df['sample_type'].astype(str)

    # Build the binary matrix: mark as presence (1) if support_percent > 0
    df['presence'] = df['support_percent'].apply(lambda x: 1 if x > 0 else 0)

    # Create the pivot table (binary matrix)
    binary_matrix = df.pivot_table(index='isoform_id',
                                   columns='protocol',
                                   values='presence',
                                   aggfunc='max',
                                   fill_value=0)
    
    # Save the binary matrix
    binary_file_path = os.path.join(binary_output_folder, os.path.basename(file_name).replace(".tsv", "_binary_matrix.txt"))
    binary_matrix.to_csv(binary_file_path, sep="\t")
    print(f"📝 Saved binary matrix: {binary_file_path}")

    # Compute Jaccard similarity matrix
    def jaccard_similarity(col1, col2):
        intersection = (col1 & col2).sum()
        union = (col1 | col2).sum()
        return intersection / union if union != 0 else 0

    protocols = binary_matrix.columns
    jaccard_df = pd.DataFrame(index=protocols, columns=protocols)

    for p1 in protocols:
        for p2 in protocols:
            jaccard_df.loc[p1, p2] = jaccard_similarity(binary_matrix[p1], binary_matrix[p2])

    jaccard_df = jaccard_df.astype(float)

    # 🔺 Plot lower triangle Jaccard heatmap to avoid repetition
    mask = np.triu(np.ones_like(jaccard_df, dtype=bool), k=1)

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        jaccard_df,
        annot=True,
        mask=mask,
        cmap="RdBu_r",           # Blue → White → Red
        vmin=0,
        vmax=1,
        fmt=".2f",
        square=True,
        linewidths=0.5,
        linecolor='white',
        cbar_kws={"label": "Jaccard Similarity"}
    )
    plt.title(f"Jaccard Similarity: {os.path.basename(file_name)}")
    plt.tight_layout()
    output_path = os.path.join(output_folder, os.path.basename(file_name).replace(".tsv", "_jaccard.png"))
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"✅ Saved Jaccard heatmap: {output_path}")

# Process each TSV file
for file in tsv_files:
    try:
        df = pd.read_csv(file, sep="\t")
        required_cols = {'Sample', 'Isoform', 'Gene', 'protocol', 'support_percent', 'sample_type'}
        if required_cols.issubset(df.columns):
            compute_jaccard(df, file)
        else:
            print(f"⚠️ Skipped {file} – missing required columns.")
    except Exception as e:
        print(f"❌ Error processing {file}: {e}")
