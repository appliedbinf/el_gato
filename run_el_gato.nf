#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads_dir = false
params.threads = 1
params.out = 'el_gato_out'

process RUN_EL_GATO_READS {
  
  conda "-c bioconda -c appliedbinf elgato=0.1.0"
  cpus 1

  input:
    tuple val(sampleId), file(reads)

  output:
    path '*_mlst.txt', emit: files

  script:
  
  r1 = reads[0]
  r2 = reads[1]

  """
  el_gato.py \
  -1 $r1 \
  -2 $r2 \
  -o out \
  -t ${task.cpus} \
  -w > ${sampleId}_mlst.txt
  """
}

process CAT {
  publishDir params.out, mode: 'copy', overwrite: true
  input:
    path files

  output:
    path 'all_mlst.txt'
  
  """
  printf "Sample\tST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA_neuAH\tmompS_min_coverage\tmompS_reads_with_primer\n" > all_mlst.txt
  cat $files >> all_mlst.txt
  """
}

workflow {
  if (params.reads_dir) {

  readPairs = Channel.fromFilePairs(params.reads_dir + "/*R{1,2}*.fastq.gz", checkIfExists: true)

  files = RUN_EL_GATO_READS(readPairs).collect()
  CAT(files)

  } else {
    print "Please provide the path to a directory containing paired reads using --reads_dir."
  }


}
