#!/usr/bin/env Rscript
#Ponsero Alise 2024
#Version  DRAFT
#This script takes a Kraken report file or a directory containing Kraken report 
# files and wrangle them into an easy to parse data format.
########## 
#usage: kraken2table.R [-h] [-i INPUT] [-t INPUT_TYPE] [-n FILES_NAMES] [-p]
#                    [--verbose] [--quietly]
#Kraken2 reports to table parser
#required arguments:
#  -i, --input : Kraken2 report input file
#optional arguments:
#  -t,  --input_type : Type of input profided: file|dir [default=file]
#  -n,  --files_names : Extensions of the Kraken files to parse [default=_profiles.txt]
#  -p,  --pct : Use the percentages instead of read counts [default=FALSE]
#other arguments:
#  -h, --help : show this help message and exit
#  -v,  --verbose : Verbose output [default=TRUE]
#  -q,  --quietly : Quiet output [default=FALSE]
##########

##########
# Libraries
suppressPackageStartupMessages(library("argparse"))
##########


##########
# create parser object
parser <- ArgumentParser()
# Required arguments
parser$add_argument("-i", "--input", type="character", required=TRUE, 
                    help="path to Kraken2 report file or to directory to parse",
                    metavar="input")
# Optional arguments
parser$add_argument("-t", "--input_type", type="character", default="file", 
                    help="Type of input profided: file|dir [default=file]",
                    metavar="input_type")
parser$add_argument("-n", "--files_names", type="character", default="_profiles.txt", 
                    help="extensions of the Kraken files to parse [default=_profiles.txt]",
                    metavar="files_names")
parser$add_argument("-p", "--pct", action="store_true", default=FALSE,
                    help="Use the percentages instead of read counts [default=FALSE]")
# Other arguments
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")

args <- parser$parse_args()
taxlevels="D,P,C,O,F,G,S"
##########

##########
# Input values sanity check
if ( ! args$input_type %in% c("dir","file") ) {
  write("INPUT ERROR: --input_type requires either file or dir as input value", stderr())
  quit(status=1)
}

if ( args$input_type == "file" ) {
  if ( ! file.exists(args$input) ) {
    write("INPUT ERROR: input file doesn't exists", stderr())
    quit(status=1)
  }
}

if ( args$input_type == "dir" ) {
  
  files2parse <- sort(list.files(args$input, pattern=args$files_names, full.names = TRUE))
  
  if ( length(files2parse)==0 ) {
    write("INPUT ERROR: input files not found in provided directory", stderr())
    quit(status=1)
  }
}

##########
# print some progress messages to stdout if "quietly" wasn't requested
if ( args$verbose ) { 
  write("Starting the parsing process...\n", stdout()) 
}

##########
### Single file parser
if ( args$input_type == "file" ) { 
  
  if ( args$verbose ) { 
    write(paste("Processing one unique file", args$input, sep=": "), stdout()) 
  }
  
  ## change to something smarter below!
  sampleID <- gsub(args$files_names, "", args$input)
  
  # load report
  report <- read.table(file = args$input, header = FALSE, sep = "\t")
  colnames(report) <- c("percent","reads_total","reads","rank","taxid","taxon")
  
  taxonomy_standard_labels <- unlist(strsplit(taxlevels,","))
  
  # prepare empty result data frame
  headers <- c("rank","taxid","taxon","taxonomy_rank","taxonomy_standard")
  headers <- c(headers, taxonomy_standard_labels)
  df = data.frame(matrix(nrow = 0, ncol = length(headers)))
  colnames(df) = headers
  
  # go through each line of report and aggregate information
  taxonomy_rank <- c("root")
  for (i in 1:nrow(report)) {
    tax <- report$taxon[i]
    taxnospaces <- gsub("^ *", "", tax)
    
    # (1) extract full taxonomy with taxa ranks
    nspaces = nchar(tax)- nchar(taxnospaces)
    indent = nspaces/2 +1 #+1 is to catch also entries without any indent, e.g. "root"
    # if this is a lower taxonomic rank, add it, otherwise reset and add
    if ( indent > length(taxonomy_rank) ) {
      taxonomy_rank <- c(taxonomy_rank,paste0(report$rank[i],"__",taxnospaces))
    } else {
      taxonomy_rank <- taxonomy_rank[1:indent-1]
      taxonomy_rank <- c(taxonomy_rank,paste0(report$rank[i],"__",taxnospaces))
    }
    taxonomy_rank_string <- paste(taxonomy_rank,collapse=";")
    
    # (2) filter taxonomy_rank to only contain entries with taxonomy_standard_labels+"__"
    taxonomy_standard = list()
    taxonomy_ranks = gsub( "__.*", "", taxonomy_rank )
    for (x in taxonomy_standard_labels) {
      taxonomy_ranks_match <- match(x, taxonomy_ranks)
      if ( !is.na(taxonomy_ranks_match) ) {
        taxonomy_clean <- gsub( ".*__", "", taxonomy_rank[taxonomy_ranks_match] )
        taxonomy_standard <- c(taxonomy_standard,taxonomy_clean)
      } else {
        taxonomy_standard <- c(taxonomy_standard,"")
      }
    }
    taxonomy_standard_string <- paste(taxonomy_standard,collapse=";")
    names(taxonomy_standard) <- taxonomy_standard_labels
    
    # (3) propagate all in results data frame
    results <- c(
      rank=report$rank[i],
      taxid=report$taxid[i],
      taxon=taxnospaces,
      taxonomy_rank=taxonomy_rank_string,
      taxonomy_standard=taxonomy_standard_string,
      taxonomy_standard)
    df <- rbind(df, results)
  }
  
  # merge with reads or pcts
  if (!args$pct) {
    counts_sel <- subset(report, select = c(taxid,reads))
  } else {
    counts_sel <- subset(report, select = c(taxid,percent))
  }
  df_counts <- merge(counts_sel, df, by.x="taxid", by.y="taxid", all.x=TRUE, all.y=TRUE)
  
  write.table(df_counts, file = paste(sampleID, ".parsed.tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")

  if ( args$verbose ) { 
    write(paste("Output writen", paste(sampleID, ".parsed.tsv", sep=""), sep=": "), stdout()) 
  }
}
##########

##########
### Directory file parser
if ( args$input_type == "dir" ) { 
  
  files2parse <- sort(list.files(args$input, pattern=args$files_names, full.names = TRUE))
  
  if ( args$verbose ) { 
    write(paste("Processing a directory", args$input, sep=": "), stdout()) 
    write(paste("found", as.character(length(files2parse)), "to process", sep=" "), stdout()) 
  }
  
  for(myfile in files2parse) {
    
    sampleID <- gsub(args$files_names, "", myfile)
    
    ## load report
    report <- read.table(file = myfile, header = FALSE, sep = "\t")
    colnames(report) <- c("percent","reads_total","reads","rank","taxid","taxon")
    
    taxonomy_standard_labels <- unlist(strsplit(taxlevels,","))
    
    # prepare empty result data frame
    headers <- c("rank","taxid","taxon","taxonomy_rank","taxonomy_standard")
    headers <- c(headers, taxonomy_standard_labels)
    df = data.frame(matrix(nrow = 0, ncol = length(headers)))
    colnames(df) = headers
    
    # go through each line of report and aggregate information
    taxonomy_rank <- c("root")
    for (i in 1:nrow(report)) {
      tax <- report$taxon[i]
      taxnospaces <- gsub("^ *", "", tax)
      
      # (1) extract full taxonomy with taxa ranks
      nspaces = nchar(tax)- nchar(taxnospaces)
      indent = nspaces/2 +1 #+1 is to catch also entries without any indent, e.g. "root"
      # if this is a lower taxonomic rank, add it, otherwise reset and add
      if ( indent > length(taxonomy_rank) ) {
        taxonomy_rank <- c(taxonomy_rank,paste0(report$rank[i],"__",taxnospaces))
      } else {
        taxonomy_rank <- taxonomy_rank[1:indent-1]
        taxonomy_rank <- c(taxonomy_rank,paste0(report$rank[i],"__",taxnospaces))
      }
      taxonomy_rank_string <- paste(taxonomy_rank,collapse=";")
      
      # (2) filter taxonomy_rank to only contain entries with taxonomy_standard_labels+"__"
      taxonomy_standard = list()
      taxonomy_ranks = gsub( "__.*", "", taxonomy_rank )
      for (x in taxonomy_standard_labels) {
        taxonomy_ranks_match <- match(x, taxonomy_ranks)
        if ( !is.na(taxonomy_ranks_match) ) {
          taxonomy_clean <- gsub( ".*__", "", taxonomy_rank[taxonomy_ranks_match] )
          taxonomy_standard <- c(taxonomy_standard,taxonomy_clean)
        } else {
          taxonomy_standard <- c(taxonomy_standard,"")
        }
      }
      taxonomy_standard_string <- paste(taxonomy_standard,collapse=";")
      names(taxonomy_standard) <- taxonomy_standard_labels
      
      # (3) propagate all in results data frame
      results <- c(
        rank=report$rank[i],
        taxid=report$taxid[i],
        taxon=taxnospaces,
        taxonomy_rank=taxonomy_rank_string,
        taxonomy_standard=taxonomy_standard_string,
        taxonomy_standard)
      df <- rbind(df, results)
    }
    
    # merge with reads or pct
    if (!args$pct) {
      counts_sel <- subset(report, select = c(taxid,reads))
    } else {
      counts_sel <- subset(report, select = c(taxid,percent))
    }
    df_counts <- merge(counts_sel, df, by.x="taxid", by.y="taxid", all.x=TRUE, all.y=TRUE)
    
    write.table(df_counts, file = paste(sampleID, ".parsed.tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
    
  }
  
  if ( args$verbose ) { 
    write(paste("Outputs writen in", args$input, sep=" "), stdout()) 
  }
  
}
