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
    path 'inputs.json', emit: f

  """
  #!/usr/bin/env python

  import sys
  import json

  inputs = {
    'read1' : "${params.read1}",
    'read2' : "${params.read2}",
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
    'logging_buffer_message' : "${params.logging_buffer_message}",
    'spades' : "${params.spades}"
  }

  # Correct input strings for not-provided inputs
  for k, v in inputs.items():
    if v == "False":
      inputs[k] = False

  with open("inputs.json", "w") as fout:
    json.dump(inputs, fout, ensure_ascii=False, indent=4)

  """

}

process UPDATE_DB {
  input:
    path f
    val db

  output:
    path f, emit: f

  """
  #!/usr/bin/env python

  import sys
  import json

  with open("${f}") as fin:
    inputs = json.load(fin)
  old_db_path = inputs['sbt']
  inputs['sbt'] = "$db"
  inputs['profile'] = inputs['profile'].replace(old_db_path, "$db")
  with open("${f}", "w") as fout:
    json.dump(inputs, fout, ensure_ascii=False, indent=4)

  """
}

process UPDATE_ASSEMBLY {
  input:
    path f
    val assembly

  output:
    path f, emit: f

  """
  #!/usr/bin/env python

  import sys
  import json

  with open("${f}") as fin:
    inputs = json.load(fin)
  inputs['assembly'] = "$assembly"
  with open("${f}", "w") as fout:
    json.dump(inputs, fout, ensure_ascii=False, indent=4)

  """
}

process UPDATE_READS {
  input:
    path f
    val read1
    val read2

  output:
    path f, emit: f

  """
  #!/usr/bin/env python

  import sys
  import json

  with open("${f}") as fin:
    inputs = json.load(fin)
  inputs['read1'] = "$read1"
  inputs['read2'] = "$read2"
  with open("${f}", "w") as fout:
    json.dump(inputs, fout, ensure_ascii=False, indent=4)

  """
}

process READS_AND_ASSEMBLY_PIPELINE {
  input:
    path f

  output:
    stdout

  """
  #!/usr/bin/env python
  
  import json
  import os

  from el_gato import el_gato

  with open("$f") as fin:
    inputs = json.load(fin)

  # Use el_gato function to abort if files missing
  el_gato.check_files(inputs)

  # Change ref file location to DB location
  ref = el_gato.Ref
  ref.file = os.path.join(inputs['sbt'], ref.file)

  # Build blast databases and other files
  el_gato.validate_ref(inputs, ref)

  # Run read + assembly pipeline
  mompS_allele = el_gato.call_momps_mapping(inputs, r1=inputs["read1"], r2=inputs["read2"],
                                    threads=inputs['threads'], ref_file=ref.file,
                                    outfile=os.path.join(inputs["sbt"], inputs["sample_name"]))
  if mompS_allele == "-":
      mompS_allele = el_gato.call_momps_pcr(inputs, assembly_file=inputs["assembly"],
                                    db=os.path.join(inputs["sbt"], "mompS" + inputs["suffix"]))
  alleles = el_gato.blast_non_momps(inputs, assembly_file=inputs["assembly"])
  alleles["mompS"] = mompS_allele

  print(el_gato.print_table(inputs,  ref, alleles))

  """
}

process ASSEMBLY_ONLY_PIPELINE {
  input:
    path f

  output:
    stdout

  """
  #!/usr/bin/env python
  
  import json
  import os

  from el_gato import el_gato

  with open("$f") as fin:
    inputs = json.load(fin)

  # Use el_gato function to abort if files missing
  el_gato.check_files(inputs)

  # Change ref file location to DB location
  ref = el_gato.Ref
  ref.file = os.path.join(inputs['sbt'], ref.file)

  # Build blast databases and other files
  el_gato.validate_ref(inputs, ref)

  # Run assembly-only pipeline
  alleles = el_gato.blast_non_momps(inputs, assembly_file=inputs["assembly"])
  alleles["mompS"] = el_gato.call_momps_pcr(inputs, assembly_file=inputs["assembly"],
                                            db=os.path.join(
                                                            inputs["sbt"], 
                                                            "mompS" + inputs["suffix"]
                                                            )
                                            )

  print(el_gato.print_table(inputs,  ref, alleles))

  """
}

process READS_ONLY_PIPELINE {
  input:
    path f

  output:
    stdout

  """
  #!/usr/bin/env python
  
  import json
  import os

  from el_gato import el_gato
  from el_gato.el_gato import SAM_data

  with open("$f") as fin:
    inputs = json.load(fin)


  # Use el_gato function to abort if files missing
  el_gato.check_files(inputs)

  # Change ref file location to DB location
  ref = el_gato.Ref
  ref.file = os.path.join(inputs['sbt'], ref.file)

  # Build blast databases and other files
  el_gato.validate_ref(inputs, ref)

  # Choose reads only or assembly followed by reads analysis
  if inputs['spades']:
    # Run read + assembly pipeline

    inputs = el_gato.genome_assembly(inputs, r1=inputs["read1"], r2=inputs["read2"],
                          out=os.path.join(inputs["sbt"], "run_spades"))

    mompS_allele = el_gato.call_momps_mapping(inputs, r1=inputs["read1"], r2=inputs["read2"],
                                      threads=inputs['threads'], ref_file=ref.file,
                                      outfile=os.path.join(inputs["sbt"], inputs["sample_name"]))
    if mompS_allele == "-":
        mompS_allele = el_gato.call_momps_pcr(inputs, assembly_file=inputs["assembly"],
                                      db=os.path.join(inputs["sbt"], "mompS" + inputs["suffix"]))
    alleles = el_gato.blast_non_momps(inputs, assembly_file=inputs["assembly"])
    alleles["mompS"] = mompS_allele

  else:
    # Run reads only pipeline
    alleles = el_gato.run_stringmlst(
      r1=inputs["read1"],
      r2=inputs["read2"],
      sbt=inputs["sbt"]
      )
    alleles["mompS"] = el_gato.check_mompS_alleles(
      r1=inputs["read1"],
      r2=inputs["read2"],
      outdir=inputs['out_prefix'],
      threads=inputs['threads'],
      ref=ref,
      db=inputs['sbt']
      )
    
  print(el_gato.print_table(inputs,  ref, alleles))

  """
}



workflow reads_and_assembly {
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
  // new_paths.db.view { "$it is in db" }
  // new_paths.assembly.view { "$it is in assembly" }
  // new_paths.r1.view { "$it is in r1" }
  // new_paths.r2.view { "$it is in r2" }

  // unpack paths to separate variables
  db = new_paths.db.view {}
  assembly = new_paths.assembly.view {}
  r1 = new_paths.r1.view {}
  r2 = new_paths.r2.view {}

  // Update inputs object with path of copied files
  f = READ_INPUTS()
  f = UPDATE_DB(f, db)
  f = UPDATE_ASSEMBLY(f, assembly)
  f = UPDATE_READS(f, r1, r2)

  READS_AND_ASSEMBLY_PIPELINE(f) | view

}


workflow assembly_only {

  to_copy = channel.of(params.sbt, params.assembly)
  destinations = channel.of("db", "assembly.fna")
  
  cp_out = CP_THING(to_copy, destinations)
  cp_out
    .branch {
        db: it =~ /db$/
        assembly: it =~ /assembly.fna$/
    }
    .set { new_paths }

  db = new_paths.db.view {}
  assembly = new_paths.assembly.view {}

  f = READ_INPUTS()
  f = UPDATE_DB(f, db)
  f = UPDATE_ASSEMBLY(f, assembly)

  ASSEMBLY_ONLY_PIPELINE(f) | view



}


workflow reads_only {
  to_copy = channel.of(params.sbt, params.read1, params.read2)
  destinations = channel.of("db", "read1.fastq.gz", "read2.fastq.gz")

  cp_out = CP_THING(to_copy, destinations)
  cp_out
    .branch {
        db: it =~ /db$/
        r1: it =~ /read1.fastq.gz$/
        r2: it =~ /read2.fastq.gz$/
    }
    .set { new_paths }

  db = new_paths.db.view {}
  r1 = new_paths.r1.view {}
  r2 = new_paths.r2.view {}

  f = READ_INPUTS()
  f = UPDATE_DB(f, db)
  f = UPDATE_READS(f, r1, r2)
  READS_ONLY_PIPELINE(f) | view
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
