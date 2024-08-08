# Metaphlan parser

## Dependencies:
- Python 3.8+
- pandas
- numpy

## How to run:
Open the Python file and change the variables at top of scirpt for your specific needs.
```
METAPHLAN_DIRECTORY="directory containing the outputs from metaphlan"
TAXA_LEVEL="The desired level you wish to concat the reads"
REL_ABUN_THRESHOLD="Threshold to remove results (by default it does not filter. 10**-3 is recommended if you want to include filtering)."
OUTFILE="Path and name for output file"
```
Then run file:
```
python read_metaphlan.py
```