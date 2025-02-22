nextflow.enable.dsl = 2

workflow {
  mutations_pre = preprocess_mutations(
    Channel.fromPath("bin/preprocess_mutations.py"),
    Channel.fromPath("data/final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf")
  )

  matrix = generate_matrix(
    Channel.fromPath("bin/generate_matrix.py"),
    mutations_pre,
  )

  extract_signatures(
    Channel.fromPath("bin/extract_signatures.py"),
    matrix
  )
}

process preprocess_mutations {
  input:
  path script
  path mutations_raw

  output:
  path "matrix"

  script:
  """
  mkdir -p ${projectDir}/results
  mkdir -p \$PWD/matrix
  python ${script} ${mutations_raw} matrix/preprocessed_mutations.maf ${projectDir}/results
  """
}

process generate_matrix {
  input:
  path script
  path preprocessed_data

  output:
  path "matrix/output/SBS/test.SBS96.all"

  script:
  """
  python ${script} ${preprocessed_data}
  """
}

process extract_signatures {
  input:
  path script
  path matrix

  output:
  path "signature_plots"

  script:
  """
  python ${script} ${matrix} signature_plots
  """
}
