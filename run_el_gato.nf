#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads_dir = false
params.assembly_dir = false
params.threads = 1
params.depth = 10
params.length = 0.3
params.sequence = 95.0
params.out = 'el_gato_out'
params.samfile = false
params.help = false

if ( params.help ) {
    help = """run_el_gato.nf: A nextflow pipeline that runs el_gato.py
             |Required arguments:
             |  Please specify either reads_dir OR assembly_dir
             |  --reads_dir      Location of the directory containing Illumina
             |                   paired-end samples. Sequence files can be gzipped.
             |  --assembly_dir   Location of the directory contaning Legionella pneumophila
             |                   genome assemblies. 
             |Optional arguments:
             |  --help           Show this help message and exit
             |  --threads        Number of threads per process (e.g. that el_gato will use per sample)
             |                   [default: ${params.threads}]
             |  --depth          Specify the minimum depth used to identify loci in paired-end reads 
             |                   [default: ${params.depth}]
             |  --length         Specify the BLAST hit length threshold for identifying multiple loci 
             |                   in assembly [default: ${params.length}]
             |  --sequence       Specify the BLAST hit percent identity threshold for identifying 
             |                   multiple loci in assembly [default: ${params.sequence}]
             |  --out            Output folder name [default: ${params.out}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

process RUN_EL_GATO_READS {
  
  cpus params.threads
  publishDir params.out, mode: 'copy', overwrite: true, pattern: '*_out/*'
  label 'elgato'

  input:
  tuple val(sampleId), file(reads)

  output:
  path '*_out/*', emit: files


  script:
  
  r1 = reads[0]
  r2 = reads[1]
  if (params.samfile == true) {
    
  """
  mkdir ${sampleId}_out/

  el_gato.py \
  -1 $r1 \
  -2 $r2 \
  -o ${sampleId}_out \
  -t ${task.cpus} \
  -d $params.depth \
  -m \
  -w > mlst.txt

  mv mlst.txt ${sampleId}_out/

  for file in \$(ls ${sampleId}_out/); do
  mv ${sampleId}_out/\$file ${sampleId}_out/${sampleId}_\$file
  done

  """
}else{
  """
  mkdir ${sampleId}_out/

  el_gato.py \
  -1 $r1 \
  -2 $r2 \
  -o ${sampleId}_out \
  -t ${task.cpus} \
  -d $params.depth \
  -w > mlst.txt

  mv mlst.txt ${sampleId}_out/

  for file in \$(ls ${sampleId}_out/); do
  mv ${sampleId}_out/\$file ${sampleId}_out/${sampleId}_\$file
  done

  """
  }
}

process RUN_EL_GATO_ASSEMBLIES {
  cpus 1
  publishDir params.out, mode: 'copy', overwrite: true, pattern: '*_out/*'
  label 'elgato'

  input:
  path assembly
  float length
  float sequence

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
  -l $params.length \
  -q $params.sequence \
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

process FINAL_JSON {
  publishDir params.out, mode: 'copy', overwrite: true
  input:
  path files

  output:
  path 'report.json'
  
  """
  echo "[" > report.tmp
  cat \$(ls *.json | head -1) >> report.tmp
  ls *.json | tail -n+2 | while read jfile; do
  echo "," >> report.tmp;
  cat \$jfile >> report.tmp;
  done
  echo "]" >> report.tmp
  mv report.tmp report.json
  """
}

process FINAL_REPORT {
  publishDir params.out, mode: 'copy', overwrite: true
  input:
  path files
  
  output:
  path 'report.pdf'
    
  """
  elgato_report.py -i *.json -o report.pdf
  """
}

workflow {
  if (params.reads_dir) {

    readPairs = Channel.fromFilePairs(params.reads_dir + "/*R{1,2}*.fastq*", checkIfExists: true)

    files = RUN_EL_GATO_READS(readPairs).collect()
    CAT(files)
    FINAL_JSON(files)
    FINAL_REPORT(files)


  } else {
    if (params.assembly_dir) {

      assemblies = Channel.fromPath(params.assembly_dir + '/*', checkIfExists: true)

      files = RUN_EL_GATO_ASSEMBLIES(assemblies).collect()
      CAT(files)
      FINAL_JSON(files)
      FINAL_REPORT(files)

    } else {
      print "Please provide the path to a directory containing paired reads using --reads_dir or the path to a directory containing assemblies using --assembly_dir."
    }
  }
}
