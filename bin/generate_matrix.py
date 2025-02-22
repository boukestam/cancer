import sys
import os
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# Read input arguments from Nextflow
input_dir = sys.argv[1]

genInstall.install("GRCh37")

matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh37", input_dir, plot=False)

print("Matrix generated")