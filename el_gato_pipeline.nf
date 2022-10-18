#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CHECK_INPUTS {

  output:
    stdout

  """
  #!/usr/bin/env python
  from el_gato import el_gato
  """

}



workflow {
  CHECK_INPUTS | view
}
