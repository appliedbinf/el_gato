#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process READ_INPUTS {

  output:
    path 'inputs.json'

  """
  #!/usr/bin/env python

  import json

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

  with open("inputs.json", "w") as fout:
    json.dump(inputs, fout, ensure_ascii=False, indent=4)

  """

}

process VALIDATE_REF {
  input:
    path 'inputs.json'

  output:
    stdout

  // """
  // cat 'inputs.json'
  // """

  """
  #!/usr/bin/env python
  
  import json

  with open('inputs.json') as fin:
    inputs = json.load(fin)

  print("printing from VALIDATE_REF")
  print(inputs)
  """


}



workflow {
  inputs = READ_INPUTS()
  VALIDATE_REF(inputs) | view
}
