import sys
import os
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
import torch

# Read input arguments from Nextflow
input_file = sys.argv[1]
output_dir = sys.argv[2]

if __name__ == "__main__":
  # Ensure output directory exists
  os.makedirs(output_dir, exist_ok=True)

  torch.set_default_tensor_type(torch.cuda.FloatTensor)
  print("PyTorch is using:", "GPU" if torch.cuda.is_available() else "CPU")

  # Run SigProfilerExtractor to find mutation signatures
  sig.sigProfilerExtractor(
    input_type="matrix",
    output=output_dir,
    input_data=input_file,
    context_type="96",
    minimum_signatures=1,
    maximum_signatures=20,
    nmf_replicates=100,
    min_nmf_iterations=10000,
    max_nmf_iterations=200000,
    matrix_normalization="gmm",
    gpu=True,
    batch_size=10
  )

  print(f"Mutation signature extraction complete. Results saved in {output_dir}")