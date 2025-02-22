import pandas as pd
import sys
import os
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import matplotlib.pyplot as plt

# Read input arguments from Nextflow
input_file = sys.argv[1]
output_file = sys.argv[2]
result_dir = sys.argv[3]

# Keep only needed columns
correct_columns = [
  'Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position', 'Strand',
  'Variant_Classification', 'Variant_Type', 'Reference_Allele',
  'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
  'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
  'Project_Code', 'Donor_ID'
]

dtype_dict = {
  'Hugo_Symbol': 'category',  # Gene symbols are repeated, so 'category' saves space
  'Chromosome': 'category',  # Few unique values (e.g., '1'-'22', 'X', 'Y')
  'Start_position': 'int32',  # Genomic coordinates are large numbers but fit in int32
  'End_position': 'int32',    # Similar to Start_position
  'Strand': 'category',  # Usually '+' or '-', so 'category' is efficient
  'Variant_Classification': 'category',  # Limited unique values, good for 'category'
  'Variant_Type': 'category',  # Few unique values, so 'category' saves memory
  'Reference_Allele': 'category',  # Typically 'A', 'T', 'C', 'G', or indels
  'Tumor_Seq_Allele1': 'category',  # Same reasoning as Reference_Allele
  'Tumor_Seq_Allele2': 'category',  # Same reasoning as Reference_Allele
  'dbSNP_RS': 'string',  # dbSNP IDs (e.g., rs1234) are best stored as strings
  'dbSNP_Val_Status': 'category',  # Few unique values
  'Tumor_Sample_Barcode': 'string',  # Unique identifiers should be stored as strings
  'Project_Code': 'category',  # Limited unique values, best as 'category'
  'Donor_ID': 'string'  # Unique patient/sample IDs, best stored as strings
}

# Load mutation data
df = pd.read_csv(input_file, sep="\t", comment="#", usecols=correct_columns, dtype=dtype_dict)

# Insert missing columns
df.insert(1, "Entrez_Gene_Id", "")
df.insert(2, "Center", "")
df.insert(3, "NCBI_Build", "")

# Remove synonymous mutations
df = df[df['Variant_Classification'] != 'Silent']

# Select a subset random tumor samples for testing
#selected_samples = df['Tumor_Sample_Barcode'].drop_duplicates().sample(n=250, random_state=42)
#df = df[df['Tumor_Sample_Barcode'].isin(selected_samples)]

# Count and plot the number of donors per cancer type
donors_by_cancer_counts = df.groupby("Project_Code")["Donor_ID"].nunique()

plt.figure(figsize=(12, 6))
donors_by_cancer_counts.sort_values(ascending=False).plot(kind="bar")
plt.xlabel("Cancer Type")
plt.ylabel("Number of Unique Samples")
plt.title("Number of Unique Samples per Cancer Type")
plt.xticks(rotation=45, ha="right")
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.savefig(os.path.join(result_dir, "donors_by_cancer_counts.png"), dpi=300, bbox_inches="tight")
plt.show()

# Keep only one cancer type
df = df[df["Project_Code"] == "Skin-Melanoma"]

# Save preprocessed file
df.to_csv(output_file, sep="\t", index=False)
print(f"Preprocessed data saved to {output_file}")
