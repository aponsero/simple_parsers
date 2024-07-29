# Centrifuge Report Parser

## What does this parser do?

This R script processes Centrifuge metagenomic classification reports. For each report file in a specified directory, it performs the following operations:

1. Reads the Centrifuge report file.
2. Performs taxonomy lookups using the TaxonomizR database, and provide a most up-to-date taxonomy classification of the reads.
3. Identifies and excludes unclassified entries (in particular taxID that are not included in NCBI anymore).
4. Aggregates read counts at each taxonomic levels.
5. Generates two output files for each input:
   - An "excluded_taxa.tsv" file containing entries that could not be classified in the new NCBI taxonomy. This happens as the latest Centrifuge database was based on an older NCBI taxonomy database, and some RaxID have been changed, merged or excluded. These TaxID will be aggregated as "unclassified" by this script.
   - A "parsed_centrifuge.tsv" file containing the final aggregated data.

## How do I install it? 

### Required R Packages

To use this script, you need to install the following R packages:

```r
install.packages(c("tidyverse", "taxonomizr"))
```

Note that the current script has been developped and tested using the following versions:

* R 4.3.3
* TaxonomizR 0.10.6
* tidyverse 2.0.0

### TaxonomizR Database

You also need to set up the TaxonomizR NCBI database to retrieve the most up to date taxonomic classifications. 

For QIB users, the TaxonomizR database is made available by the QIB Core Bioinformatics. The database is available at :

```
/qib/platforms/Informatics/transfer/outgoing/databases/taxonomizR
```

For other users, you want to install your own database on your system :

```r
# Create a directory for the database:
mkdir TaxonomizR_db

# Download and prepare the database (this may take some time):
library(taxonomizr)

# Set the path for your database
db_path <- "TaxonomizR_db/NCBItax.sql"

# Download and prepare the database
prepareDatabase(db_path)
```

Note: The database preparation might take several hours and require significant disk space (>100GB).

## How do I use this parser ?

If you can easily mount the database (mac systems) or have the taxonomizR database on your computer, you can download and save the script as centrifuge2aggregation.R.

In your R environment or script, source the file:

```r
library(tidyverse)
library(taxonomizr)
source("path/to/centrifuge2aggregation.R")

# Call the process_centrifuge_files function with your directory and database paths:

directory_path <- "path/to/your/centrifuge/reports"
ncbi_db_path <- "path/to/your/TaxonomizR_db/NCBItax.sql"
filepattern <- "_report.txt"
process_centrifuge_files(directory_path, ncbi_db_path, filepattern)

```
Replace "path/to/your/centrifuge/reports" with the path to the directory containing your Centrifuge report files, and "path/to/your/TaxonomizR_db/NCBItax.sql" with the path to your TaxonomizR database file.

The script will process all files ending by "report.txt" in the specified directory, creating corresponding parsed and excluded taxa files for each input file.

## Run the script on QIB cloud

To use the parser on QIB cloud, you will need to first install the tidyverse and taxonomizR packages in a conda environment: 

```bash

conda create -n taxonomizR conda-forge::r-tidyverse conda-forge::r-taxonomizr

conda activate taxonomizR

```

Once your conda environment is created and activated, you can use R and this parser. Download the parser from github:

```bash

wget https://raw.githubusercontent.com/aponsero/simple_parsers/main/Centrifuge2aggregation/centrifuge2aggregation.R

```

Then you can use the parser as above: 

```r
library(tidyverse)
library(taxonomizr)
source("path/to/centrifuge2aggregation.R")

# Call the process_centrifuge_files function with your directory and database paths:

directory_path <- "path/to/your/centrifuge/reports"
ncbi_db_path <- "path/to/your/TaxonomizR_db/NCBItax.sql"
filepattern <- "_report.txt"
process_centrifuge_files(directory_path, ncbi_db_path, filepattern)

```

Make sure that the QIB Informatics folder is mounted in your VM by following [these instructions](https://quadram-institute.atlassian.net/wiki/spaces/BSUP/pages/19497039/How+to+access+data+on+a+QIB+Cloud+VM#How-to-mount-incoming/outgoing) in order to have access to the taxonomizR database.


