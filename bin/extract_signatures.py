import sys
import os
import pandas as pd
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig

# Read input arguments from Nextflow
input_dir = sys.argv[1]  # Path to preprocessed mutations file
output_dir = sys.argv[2]  # Path to output directory

if __name__ == "__main__":
  genInstall.install("GRCh37")

  matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh37", input_dir, plot=False)

  # Ensure output directory exists
  os.makedirs(output_dir, exist_ok=True)

  # Run SigProfilerExtractor to find mutation signatures
  sig.sigProfilerExtractor(
    input_type="matrix",
    output=output_dir,
    input_data="matrix_input/output/SBS/test.SBS96.all",
    context_type="96",
    minimum_signatures=1,
    maximum_signatures=6,
    nmf_replicates=100,
    min_nmf_iterations=10000,
    max_nmf_iterations=200000,
    matrix_normalization="gmm",
    cpu=4  # Limit CPU usage for efficiency
  )

  print(f"Mutation signature extraction complete. Results saved in {output_dir}")