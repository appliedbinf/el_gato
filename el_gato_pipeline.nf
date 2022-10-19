#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CP_THING {
  input:
    val infile
    val outfile

  output:
    path outfile

  """
  cp -r ${infile} ${outfile}
  """
}


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

  """
  #!/usr/bin/env python
  
  import json

  from el_gato import el_gato

  with open('inputs.json') as fin:
    inputs = json.load(fin)

  ref = el_gato.Ref

  print("running validate_ref()")

  el_gato.validate_ref(inputs, ref)

  """


}

process HEAD_TEST {
  input:
    path file

  output:
    stdout

  """
  head -1 ${file}
  """
}

process LS_TEST {
  input:
    path dir

  output:
    stdout

  """
  ls ${dir}
  """
}

workflow {

  if ( params.read1 != "False" && params.read2 != "False") {
    if ( params.assembly != "False" ){
      // Make list of source and destination for inputs
      to_copy = channel.of(params.sbt, params.read1, params.read2, params.assembly)
      destinations = channel.of("db", "read1.fastq", "read2.fastq", "assembly.fna")
      
      // Copy inputs to Nextflow work dir and create lookup object
      cp_out = CP_THING(to_copy, destinations)
      cp_out
        .branch {
            db: it =~ /db$/
            assembly: it =~ /assembly.fna$/
            r1: it =~ /read1.fastq$/
            r2: it =~ /read2.fastq$/
        }
        .set { new_paths }

      // Paths can now be accessed with the following
      new_paths.db.view { "$it is in db" }
      new_paths.assembly.view { "$it is in assembly" }
      new_paths.r1.view { "$it is in r1" }
      new_paths.r2.view { "$it is in r2" }
    } else {
      to_copy = channel.of(params.sbt, params.read1, params.read2)
      destinations = channel.of("db", "read1.fastq", "read2.fastq")

      cp_out = CP_THING(to_copy, destinations)
      cp_out
        .branch {
            db: it =~ /db$/
            r1: it =~ /read1.fastq$/
            r2: it =~ /read2.fastq$/
        }
        .set { new_paths }

        new_paths.db.view { "$it is in db" }
        new_paths.r1.view { "$it is in r1" }
        new_paths.r2.view { "$it is in r2" }
    }
    
  } else {
    if ( params.assembly != "False" ){

      to_copy = channel.of(params.sbt, params.assembly)
      destinations = channel.of("db", "assembly.fna")
      
      cp_out = CP_THING(to_copy, destinations)
      cp_out
        .branch {
            db: it =~ /db$/
            assembly: it =~ /assembly.fna$/
        }
        .set { new_paths }

        new_paths.db.view { "$it is in db" }
        new_paths.assembly.view { "$it is in assembly" }

    } else {
      println ( "ERROR! You must provide either reads, or an assembly, or both." )
    }
  }

  // inputs = READ_INPUTS()
  // VALIDATE_REF(inputs) | view
}
