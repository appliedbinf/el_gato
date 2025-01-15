# Reporting Module  

We provide a script that generates a PDF report of each el_gato run using the report.json file generated in the output folder for each sample.
Nextflow generates this report by default, but you must run it manually if you run el_gato for individual samples. You can generate a report for one or more samples, including assembly, and read reports.

* [Dependencies](#dependencies)
* [Usage](#usage)
  * [Quickstart](#quickstart)
  * [All available arguments](#all-available-arguments)
* [Example image of pdf report](#example-image-of-pdf-report)

## Dependencies
  * [fpdf2](https://github.com/py-pdf/fpdf2)

## Usage

### Quickstart
```
elgato_report.py -i <path/to/report1.json> [<path/to/report2.json> ...] -o <path/to/output/report.pdf>
```

### All available arguments
Usage information printed when running elgato_report.py with `-h` or `--help`

```
options:
  -h, --help            show this help message and exit
  -i, --input_jsons     path to one or more report.json files
  -o, --out_report      desired output pdf file path
  -s, --shorten_names   shorten long sample and contig names to prevent line wrapping
  -n, --no_header       Do not include the header in the report
  --custom_header       Provide a custom header as a string in your command
  --header_file         Provide a custom header in a text file
```

## Example image of pdf report
[will include when NC printing into report is resolved]



