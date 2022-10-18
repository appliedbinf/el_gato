#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CHECK_INPUTS {

  output:
    stdout

  """
  #!/usr/bin/env python
  from el_gato import el_gato

  inputs = {
    'read1' : ${params.read1},
    'read2' : ${params.read2},
    'assembly' : "${params.assembly}",
    'threads' : ${params.threads},
    'out_prefix' : "${params.out_prefix}",
    'sample_name' : "${params.sample_name}",
    'log' : "${params.log}",
    'sbt' : "${params.sbt}",
    'suffix' : "${params.suffix}",
    'profile' : "${params.profile}",
    'verbose' : ${params.verbose},
    'overwrite' : ${params.overwrite},
    'depth' : ${params.depth},
    'analysis_path' : "${params.analysis_path}",
    'logging_buffer_message' : "${params.logging_buffer_message}"
  }

  print(inputs)

  """

}



workflow {
  CHECK_INPUTS | view
}
