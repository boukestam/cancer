import pandas as pd
import sys

# Read input file (from Nextflow parameter)
input_file = sys.argv[1]
output_file = sys.argv[2]

# Load mutation data
df = pd.read_csv(input_file, sep="\t", comment="#", low_memory=False)

# Remove synonymous mutations
df = df[df['Variant_Classification'] != 'Silent']

# Select a subset random tumor samples for testing
selected_samples = df['Tumor_Sample_Barcode'].drop_duplicates().sample(n=250, random_state=42)
df = df[df['Tumor_Sample_Barcode'].isin(selected_samples)]

# Ensure correct column order
correct_columns = [
  'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome',
  'Start_Position', 'End_Position', 'Strand',
  'Variant_Classification', 'Variant_Type', 'Reference_Allele',
  'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
  'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
  'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1',
  'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
  'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1',
  'Match_Norm_Validation_Allele2', 'Verification_Status',
  'Validation_Status', 'Mutation_Status', 'Sequencing_Phase',
  'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File',
  'Sequencer', 't_ref_count', 't_alt_count', 'n_ref_count', 'n_alt_count',
  'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'RefSeq',
  'Protein_position', 'Codons', 'Hotspot'
]
df = df[correct_columns]

# Save preprocessed file
df.to_csv(output_file, sep="\t", index=False)
print(f"Preprocessed data saved to {output_file}")
