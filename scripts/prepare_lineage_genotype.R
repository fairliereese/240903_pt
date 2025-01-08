#!/usr/bin/env Rscript

#################
###  IMPORTS  ###
#################

# Import required packages
if (
  suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE
) {
  stop("[ERROR] Package 'optparse' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("dplyr"))) == FALSE
) {
  stop("[ERROR] Package 'dplyr' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("data.table"))) == FALSE
) {
  stop("[ERROR] Package 'data.table' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("tibble"))) == FALSE
) {
  stop("[ERROR] Package 'tibble' required! Aborted.")
}


#######################
###  PARSE OPTIONS  ###
#######################

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly = FALSE)[4]))

# Build description message
description <- "Prepare lineage genotype for sQTLseeker2-nf.\n"
author <- "Author: Iris Mestres"
version <- "Version: 1.1.0 (SEP-2024)"
requirements <- "Requires: optparse, dplyr, data.table, tibble"
msg <- paste(description, author, version, requirements, sep = "\n")

# Define list of arguments
option_list <- list(
  make_option(
    "--genotype",
    action = "store",
    type = "character",
    help = "Path to the compressed genotype file (.vcf.gz). Required!",
    metavar = "file"
  ),
  make_option(
    "--metadata",
    action = "store",
    type = "character",
    help = "Path to metadata file (.tsv). Required!",
    metavar = "file"
  ),
  make_option(
    "--in-012",
    action = "store",
    type = "character",
    help = "Path to the phenotype matrix in '012' format file. Required!",
    metavar = "file"
  ),
  make_option(
    "--pos-012",
    action = "store",
    type = "character",
    help = "Path to the position matrix in '012' format file. Required!",
    metavar = "file"
  ),
  make_option(
    "--indv-012",
    action = "store",
    type = "character",
    help = "Path to the individiual matrix in '012' format file. Required!",
    metavar = "file"
  ),
  make_option(
    "--out-file",
    action = "store",
    type = "character",
    default = file.path( getwd(), "lineage_in_data/" ),
    help = "Absolute filename to where output files shall be written. Default: $PWD/lineage_in_data",
    metavar = "file"
  ),
  make_option(
    c( "-h", "--help" ),
    action = "store_true",
    default = FALSE,
    help = "Show this information and die."
  ),
  make_option(
    c( "-u", "--usage" ),
    action = "store_true",
    default = FALSE,
    dest = "help",
    help = "Show this information and die."
  ),
  make_option(
    c( "-v", "--verbose" ),
    action = "store_true",
    default = FALSE,
    help = "Print log messages to STDOUT."
  )
)

# Parse command-line arguments
opt_parser <-
  OptionParser(
    usage = paste(
      "Usage:",
      script,
      "--genotype <path/to/genotype.vcf.gz>",
      "--metadata <path/to/metadata.tsv>",
      "--in-012 <path/to/matrix.012>",
      "--pos-012 <path/to/position.012.pos>",
      "--indv-012 <path/to/donor-id.012.indv>",
      "[OPTIONS]\n",
      sep = " "
    ),
    option_list = option_list,
    add_help_option = FALSE,
    description = msg
  )
opt <- parse_args( opt_parser )

# Re-assign variables
genotype.vcf.gz <- opt$`genotype`
metadata <- opt$`metadata`
in.012 <- opt$`in-012`
pos.012 <- opt$`pos-012`
indv.012  <- opt$`indv-012`
out.file <- opt$`out-file`
verb <- opt$`verbose`

# Validate required arguments
if ( is.null( genotype.vcf.gz ) | is.null( metadata) | is.null( in.012 ) | is.null( pos.012 ) | is.null( indv.012 ) ) {
  print_help(opt_parser)
  stop("[ERROR] Required argument missing! Aborted.")
}


######################
###    FUNCTIONS   ###
######################

## ASSUMPTIONS
### The input VCF file is '.gz' compressed
### The output path has the extension '.gz' so the file is compressed
### Donor IDs are just the Donor number
get_genotype_info <- function( in_vcf, pos.012 ) {
  i_tbl <- data.table::fread( in_vcf, sep = "\t",
                              header = TRUE, skip = "CHROM",
                              select = c("#CHROM", "POS", "ID", "REF", "ALT") )

  colnames( i_tbl ) <- c("chr", "pos", "snpId", 'ref', 'alt')
  i_tbl <- i_tbl %>%
    mutate(new_snpid = paste(chr, pos, ref, alt, sep = ":"))

  i_tbl$chr <- unlist( lapply( i_tbl$chr,
                               function(x) strsplit( x, "chr" )[[1]][2] ))


  i_pos.012 <- read.table( pos.012, sep="\t", col.names = c("chr", "pos"))
  i_pos.012$chr <- unlist( lapply( i_pos.012$chr,
                               function(x) strsplit( x, "chr" )[[1]][2] ))

  # TODO - remove
  write.table( i_pos.012, paste(c( 'i_pos', "temp.tsv"), collapse="_" ),
               sep = "\t", quote = FALSE,
               row.names = FALSE, col.names = TRUE )
  write.table( i_tbl.012, paste(c( 'i_tbl', "temp.tsv"), collapse="_" ),
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE )

  o_pos.012 <- i_pos.012 %>% dplyr::left_join( i_tbl, by=c( "chr", "pos" ))

# TODO -- filter on non-duplicated SNP ids

  write.table( o_pos.012, paste(c( pos.012, "extended.tsv"), collapse="_" ),
               sep = "\t", quote = FALSE,
               row.names = FALSE, col.names = TRUE )

  return ( tibble::tibble( chr=i_tbl$chr,
                           start=i_tbl$pos,
                           end=i_tbl$pos,
                           snpId=i_tbl$snpId ) )
}

get_genotype_012_format <- function( in.012 , indv.012, snpId ) {

  donors <- data.table::fread( in.012, sep = "\t", drop = 1)

  indv.012 <- unlist( lapply( indv.012$V1, function(x)
    paste( c( "popCell", x ),  collapse = "_" )))

  donors <- as.data.frame( t( donors ))
  colnames( donors ) <- indv.012

  donors$snpId <- snpId

  return( donors )
}


######################
###      MAIN      ###
######################


o.genotype.vcf.gz <- out.file

if (verb) {
  cat( "Creating genotype file...\n" )
  cat( paste( "Modifying SNPs data...\n", timestamp(), "\n" ))
}
genotype.info <- get_genotype_info( in_vcf = genotype.vcf.gz,
                                    pos.012 = pos.012 )

if (verb) {
  cat( paste( "Modifying donors data...\n", timestamp(), "\n" ))
}

indv.012 <- read.table( indv.012, sep="\t" )
indId <- data.table::fread( metadata, select = c( "indId" ))
snpId <- data.table::fread( paste(c( pos.012, "extended.tsv" ), collapse = "_" ),
                            sep = "\t", header = TRUE, select = c( "snpId" ))

genotype.donors <- get_genotype_012_format( in.012 = in.012,
                                            indv.012 = indv.012,
                                            snpId = unlist( snpId ))

genotype.donors <- genotype.donors %>% dplyr::select( "snpId", indId$indId )

if ( verb ) {
  cat( paste( "Writing genotype file ...\n", timestamp(), "\n"))
}


genotype <- dplyr::full_join(genotype.info, genotype.donors, by="snpId") %>%
  dplyr::filter(snpId != ".")

readr::write_tsv( genotype, o.genotype.vcf.gz )

if ( verb ) {
  cat( paste( "Done\n", timestamp(), "\n" ))
}
