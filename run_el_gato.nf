#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads_dir = false
params.assembly_dir = false
params.threads = 1
params.out = 'el_gato_out'

process RUN_EL_GATO_READS {
  conda "-c conda-forge -c bioconda -c appliedbinf elgato"
  cpus 1
  publishDir params.out, mode: 'copy', overwrite: true, pattern: '*_out/*'
  label 'elgato'

  input:
    tuple val(sampleId), file(reads)

  output:
    path '*_out/*', emit: files


  script:
  
  r1 = reads[0]
  r2 = reads[1]

  """
  mkdir ${sampleId}_out/

  el_gato.py \
  -1 $r1 \
  -2 $r2 \
  -o ${sampleId}_out \
  -t ${task.cpus} \
  -w > mlst.txt

  mv mlst.txt ${sampleId}_out/

  for file in \$(ls ${sampleId}_out/); do
  mv ${sampleId}_out/\$file ${sampleId}_out/${sampleId}_\$file
  done

  """
}

process RUN_EL_GATO_ASSEMBLIES {
  conda "-c conda-forge -c bioconda -c appliedbinf elgato"
  cpus 1
  publishDir params.out, mode: 'copy', overwrite: true, pattern: '*_out/*'
  label 'elgato'

  input:
    path assembly

  output:
    path '*_out/*', emit: files


  script:
  """
  sample_id=$assembly
  sample_id=\${sample_id%.*}

  el_gato.py \
  -a $assembly \
  -o \${sample_id}_out \
  -t ${task.cpus} \
  -w > mlst.txt

  mv mlst.txt \${sample_id}_out/

  for file in \$(ls \${sample_id}_out/); do
  mv \${sample_id}_out/\$file \${sample_id}_out/\${sample_id}_\$file
  done

  """
}

process CAT {
  publishDir params.out, mode: 'copy', overwrite: true
  input:
    path files

  output:
    path 'all_mlst.txt'
  
  """
  printf "Sample\tST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA_neuAH\n" > all_mlst.txt
  cat \$(ls ${files} | grep mlst | grep -v possible_mlsts) >> all_mlst.txt
  """
}

workflow {
  if (params.reads_dir) {

    readPairs = Channel.fromFilePairs(params.reads_dir + "/*R{1,2}*.fastq*", checkIfExists: true)

    files = RUN_EL_GATO_READS(readPairs).collect()
    CAT(files)


  } else {
    if (params.assembly_dir) {

      assemblies = Channel.fromPath(params.assembly_dir + '/*', checkIfExists: true)

      files = RUN_EL_GATO_ASSEMBLIES(assemblies).collect()
      CAT(files)

    } else {
      print "Please provide the path to a directory containing paired reads using --reads_dir or the path to a directory containing assemblies using --assembly_dir."
    }
  }
}
