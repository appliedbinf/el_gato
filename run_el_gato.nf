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

process RUN_EL_GATO {
  
  conda "-c bioconda -c appliedbinf elgato"
  cpus 1

  input:
    val read1
    val read2
    val assembly
    val db

  output:
    stdout

  """
  el_gato.py -h
  """
}




workflow reads_and_assembly {
  // // Make list of source and destination for inputs
  // to_copy = channel.of(params.sbt, params.read1, params.read2, params.assembly)
  // destinations = channel.of("db", "read1.fastq", "read2.fastq", "assembly.fna")
  
  // // Copy inputs to Nextflow work dir and create lookup object
  // cp_out = CP_THING(to_copy, destinations)
  // cp_out
  //   .branch {
  //       db: it =~ /db$/
  //       assembly: it =~ /assembly.fna$/
  //       r1: it =~ /read1.fastq$/
  //       r2: it =~ /read2.fastq$/
  //   }
  //   .set { new_paths }

  // // Paths can now be accessed with the following
  // // new_paths.db.view { "$it is in db" }
  // // new_paths.assembly.view { "$it is in assembly" }
  // // new_paths.r1.view { "$it is in r1" }
  // // new_paths.r2.view { "$it is in r2" }

  // // unpack paths to separate variables
  // db = new_paths.db.view {}
  // assembly = new_paths.assembly.view {}
  // r1 = new_paths.r1.view {}
  // r2 = new_paths.r2.view {}

  RUN_EL_GATO(params.read1, params.read2, params.assembly, params.sbt) | view

}

workflow assembly_only {

  // to_copy = channel.of(params.sbt, params.assembly)
  // destinations = channel.of("db", "assembly.fna")
  
  // cp_out = CP_THING(to_copy, destinations)
  // cp_out
  //   .branch {
  //       db: it =~ /db$/
  //       assembly: it =~ /assembly.fna$/
  //   }
  //   .set { new_paths }

  // db = new_paths.db.view {}
  // assembly = new_paths.assembly.view {}

  RUN_EL_GATO(params.read1, params.read2, params.assembly, params.sbt) | view
}

workflow reads_only {
  // to_copy = channel.of(params.sbt, params.read1, params.read2)
  // destinations = channel.of("db", "read1.fastq.gz", "read2.fastq.gz")

  // cp_out = CP_THING(to_copy, destinations)
  // cp_out
  //   .branch {
  //       db: it =~ /db$/
  //       r1: it =~ /read1.fastq.gz$/
  //       r2: it =~ /read2.fastq.gz$/
  //   }
  //   .set { new_paths }

  // db = new_paths.db.view {}
  // r1 = new_paths.r1.view {}
  // r2 = new_paths.r2.view {}

  RUN_EL_GATO(params.read1, params.read2, params.assembly, params.sbt) | view
}


workflow {

  if ( params.read1 != "False" && params.read2 != "False") {
    if ( params.assembly != "False" ){

      reads_and_assembly()

    } else {

      reads_only()

    }
  } else {
    if ( params.assembly != "False" ){

      assembly_only()

    } else {

      println ( "ERROR! You must provide either reads, or an assembly, or both." )

    }
  }
}
