# Reporting Module  

We provide a script that generates a PDF report for each el_gato run using the report.json file located in the output folder for each sample. You can generate reports for one or more samples, including assembly, and review the results.

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



