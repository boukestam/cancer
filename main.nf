nextflow.enable.dsl = 2

workflow {
  mutations_raw = fetch_tcga_data()

  mutations_pre = preprocess_mutations(
    Channel.fromPath("bin/preprocess_mutations.py"),
    mutations_raw,
  )

  extract_signatures(
    Channel.fromPath("bin/extract_signatures.py"),
    mutations_pre,
  )
}

process fetch_tcga_data {
  cache true

  output:
  path "mutations.tsv"

  script:
  """
  # Create a directory for data
  mkdir -p data && cd data

  # Download METABRIC dataset
  wget -q https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz

  # Extract the tar.gz file
  tar -xzf brca_metabric.tar.gz

  # Move the mutation data to the root of the data folder
  mv brca_metabric/data_mutations.txt ../mutations.tsv

  # Cleanup
  rm -rf brca_metabric brca_metabric.tar.gz
  """
}

process preprocess_mutations {
  input:
  path script
  path mutations_raw

  output:
  path "matrix_input"

  script:
  """
  mkdir -p \$PWD/matrix_input
  python ${script} ${mutations_raw} matrix_input/preprocessed_mutations.maf
  """
}

process extract_signatures {
  input:
  path script
  path preprocessed_data

  output:
  path "signature_plots"

  script:
  """
  python ${script} ${preprocessed_data} signature_plots
  """
}
