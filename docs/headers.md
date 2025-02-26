

# Headers
Headers are included in outputs for the samtools coverage command and blast results. Header definitions are as follows:

## samtools coverage headers

| Column header | Meaning 
|:-------------:|:---------------------------------------------------:|
| rname         | Locus name                                          |
| numreads      | Number reads aligned to the region (after filtering)|
| covbases      | Number of covered bases with depth >= 10             |
| coverage      | Percentage of covered bases [0..100]                |
| meandepth     | Mean depth of coverage                              |
| meanbaseq     | Mean baseQ in covered region                        |
| meanmapq      | Mean mapQ of selected reads                         | 

## BLASTn output headers

| Column header | Meaning                             |
|:-------------:|:-----------------------------------:|
| qseqid        | Query sequence id                   |
| sseqid        | Subject (matched allele) id         |
| pident        | Percentage of identical matches     |
| length        | Alignment length (sequence overlap) |
| mismatch      | Number of mismatches                |
| gapopen       | Mumber of gap openings              |
| qstart        | Start of alignment in query         |
| qend          | End of alignment in query           |
| sstart        | Start of alignment in subject       |
| send          | End of alignment in subject         |
| evalue        | Expect value                        |
| bitscore      | Bit score                           |
| qlen          | Query sequence length               |
| slen          | Subject sequence length             |
