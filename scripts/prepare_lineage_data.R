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
  suppressWarnings(suppressPackageStartupMessages(require("readr"))) == FALSE
) {
  stop("[ERROR] Package 'readr' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE
) {
  stop("[ERROR] Package 'rtracklayer' required! Aborted.")
}


#######################
###  PARSE OPTIONS  ###
#######################

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly = FALSE)[4]))

# Build description message
description <- "Prepare lineage data for sQTLseeker2-nf.\n"
author <- "Author: Iris Mestres"
version <- "Version: 1.1.0 (SEP-2024)"
requirements <- "Requires: optparse, dplyr, readr, rtracklayer"
msg <- paste(description, author, version, requirements, sep = "\n")

# Define list of arguments
option_list <- list(
  make_option(
    "--gtf",
    action = "store",
    type = "character",
    help = "Path to the annotation file (.gtf). Required!",
    metavar = "file"
  ),
  make_option(
    "--metadata",
    action = "store",
    type = "character",
    help = "Path to the metadata file (.tsv). Required!",
    metavar = "file"
  ),
  make_option(
    "--tcc-dir",
    action = "store",
    type = "character",
    help = "Path to the directory containing the TCC expression file. Required!",
    metavar = "directory"
  ),
  make_option(
    "--suffix",
    action = "store",
    type = "character",
    help = "Suffix of the TCC expression file (.rds). Required!",
    metavar = "string"
  ),
  make_option(
    "--id-dir",
    action = "store",
    type = "character",
    help = "Path to the TCC and gene IDs table (.tsv). Required!",
    metavar = "file"
  ),
  make_option(
    "--id-rds",
    action = "store",
    type = "character",
    help = "Path to the TCC and gene IDs table (.tsv). Required!",
    metavar = "file"
  ),
  make_option(
    "--pca-dir",
    action = "store",
    type = "character",
    help = "Absolute path to where PCA files shall be read from. Required!",
    metavar = "directory"
  ),
  make_option(
    "--out-dir",
    action = "store",
    type = "character",
    default = file.path( getwd(), "lineage_in_data/" ),
    help = "Absolute path to where output files shall be written. Default: $PWD/lineage_in_data",
    metavar = "directory"
  ),
  make_option(
    "--lineage",
    action = "store",
    type = "character",
    default = "MONO",
    help = "Lineages from which to create TCC expression files. Default: MONO.",
    metavar = "string"
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
      "--gtf <path/to/annotation.gtf>",
      "--metadata <path/to/metadata.tsv>",
      "--tcc-dir <path/to/tcc-expression/dir>",
      "--suffix=<suffix-of-tcc-expression.rds>",
      "--id-dir <path/to/tcc-gene-ids/dir>",
      "--id-rds=<suffix-of-tcc-gene-ids.rds>",
      "--pca-dir <path/to/pca/dir>",
      "[OPTIONS]\n",
      sep = "\n"
    ),
    option_list = option_list,
    add_help_option = FALSE,
    description = msg
  )
opt <- parse_args( opt_parser )

# Re-assign variables
annotation.gtf <- opt$`gtf`
metadata.tsv <- opt$`metadata`
tcc.exp.rds <- opt$`suffix`
tcc.exp.dir <- opt$`tcc-dir`
id.tcc.dir <- opt$`id-dir`
id.tcc.rds <- opt$`id-rds`
pca.dir <- opt$`pca-dir`
out.dir <- opt$`out-dir`
lineage <- opt$`lineage`
verb <- opt$`verbose`

# Validate required arguments
if (is.null(annotation.gtf) | is.null(metadata.tsv) |  is.null(id.tcc.dir) | is.null(id.tcc.rds) |is.null(tcc.exp.dir) | is.null(tcc.exp.rds) | is.null(pca.dir)) {
  print_help(opt_parser)
  stop("[ERROR] Required argument missing! Aborted.")
}


######################
###    FUNCTIONS   ###
######################

## ASSUMPTIONS
### Structure rds: main -> CONDITION -> TRAIT
### Gene IDs stored in 'gene_id'
### Feature IDs stored in 'feature_ids' and separated by a semi-colon ';'
get_gene_trans_ids_all <- function( in_rds ) {
  data <- readRDS( in_rds )
  conditions <- c( "COV", "IAV", "NS" )

  ids.all <- tibble( "trId" = character(),
                     "geneId" = character() )
  is.first <- TRUE

  for ( condition in conditions ) {
    data_cond <- data[[condition]][["POP"]]

    data_cond <- data_cond %>%
      dplyr::select( gene_id, feature_ids )

    ids <- tibble::tibble( "trId" = character(), "geneId" = character() )

    for ( i in seq( 1, nrow( data_cond ) )) {
      trans <- strsplit( data_cond$feature_ids[i], ';' )[[1]]

      for ( j in trans ) {
        ids <- ids %>% tibble::add_row( trId = j,
                                        geneId = data_cond$gene_id[i] )
      }
    }

    if ( is.first ) {
      ids.all <- ids
      is.first <- FALSE
    } else {
      ids.all <- ids.all %>% dplyr::full_join( ids, by = c( "geneId", "trId" ))
    }
  }
  return( ids.all )
}


## ASSUMPTIONS
### Structure rds: main -> CONDITION
### Fields in table have the format 'DONOR.ID_CONDITION'
### Rows are named by TCC IDs
### At least one TCC from 'ids' is in in_rds$CONDITION
### 'ids.all' has the fields 'trId' and 'geneId' with the TCC and gene IDs
### respectively
get_transcript_expression_file_all <- function( in_rds, ids.all, out_tbl ) {
  data <- readRDS( in_rds )
  conditions <- c( "COV", "IAV", "NS" )
  is.first <- TRUE

  for ( condition in conditions ){
    data.all <- data[[condition]]
    data.all$trId <- rownames( data.all )

    if ( is.first ) {
      o_tbl <- ids.all %>% dplyr::full_join( data.all, by="trId" )
      is.first <- FALSE
    } else {
      o_tbl <- o_tbl %>% dplyr::full_join( data.all, by="trId" )
    }
  }

  o_tbl <- o_tbl %>%
    dplyr::mutate_if( is.numeric, ~replace(., is.na(.), 0))

  readr::write_tsv( o_tbl, out_tbl )

  return( colnames( o_tbl )[3: ncol( o_tbl )] )
}

## ASSUMPTIONS
### The input annotation file is a GTF with the fields 'seqnames', 'start',
### 'end', 'type' and 'gene_id'
### 'ids' has the fields 'trId' and 'geneId' with the TCC and gene IDs
### respectively
get_gene_location_file <- function( in_gtf, ids = NULL, out_tbl ) {

  genes <- as.data.frame( rtracklayer::import( in_gtf ) )
  genes <- genes[ genes$type == "gene", ] %>%
    select( seqnames, start, end, gene_id )
  colnames( genes ) <- c("chr", "start", "end", "geneId")

  if ( !is.null( ids )){
    genes <- genes %>% dplyr::filter( geneId %in% ids$geneId )
  }
  genes$chr <- unlist(lapply(
    genes$chr, function(x) strsplit( as.character(x), "chr" )[[1]][2] ))

  write.table( genes, out_tbl,
               sep = "\t",
               col.names = TRUE,
               row.names = FALSE,
               quote = FALSE )
}


get_pc_tbl <- function( pca.path ) {
  condition <- c( "NS", "COV", "IAV" )

  for ( cond in condition ) {

    pca <- readr::read_table( pca.path, col_names = TRUE )
    pca <- pca %>%
      dplyr::select( -"#FID" ) %>%
      dplyr::mutate( IID = paste(IID, cond, sep = "_" )) %>%
      dplyr::rename( "sampleId" = IID )

    if ( cond == "NS" ) {
      pca.all <- pca
    } else {
      pca.all <- pca.all %>%
        dplyr::full_join( pca, by = c( "sampleId",
                                       "PC1" , "PC2", "PC3", "PC4", "PC5" ))
    }
  }

  return( pca.all )

}

## ASSUMPTIONS
### The fields in the file are 'Run', 'Library', 'Well', 'POP', 'Donor.ID',
### 'Condition', 'TP', 'in_final_dataset', 'ncells', 'POPASH', 'POPEUB',
### 'GenderF', 'Age', 'Mortality', 'CMV', 'Sample.ID' and 'Donor.ID_Condition'.
### At least one Donor.ID_Condition is in the transcript expression file and
### follow the same name formatting.
### The covariates to be used are Age and Run
get_metadata_file_all <- function( in.tbl, pca.tbl, out_tbl ){
  data <- read.table( in.tbl, header = TRUE, sep="\t" )

  data <- data  %>% dplyr::filter( POP != "ASH",
                                   in_final_dataset == TRUE,
                                   Run != 16 ) %>%
    dplyr::select( Donor.ID,
                   Donor.ID_Condition,
                   Condition,
                   Age, Run) %>%
    dplyr::group_by( Donor.ID_Condition ) %>%
    dplyr::arrange( Donor.ID_Condition )

  data$Donor.ID <- unlist(lapply(
    data$Donor.ID, function(x) paste( "popCell", x, sep = "_" )))

  data <- data %>% dplyr::rename( "indId" = Donor.ID,
                                  "sampleId" = Donor.ID_Condition,
                                  "group" = Condition,
                                  "cov.age" = Age,
                                  "cov.run" = Run)

  data <- data %>%
    dplyr::distinct( .keep_all = TRUE ) %>%
    dplyr::full_join( pca.tbl, by = "sampleId" )

  write.table( data, out_tbl,
               sep = "\t",
               col.names = TRUE,
               row.names = FALSE,
               quote = FALSE )
}

######################
###      MAIN      ###
######################

if (verb) {
  cat( "Setting output directory...\n" )
}

if ( !dir.exists( out.dir ) ) { dir.create( out.dir ) }

o.genes.loc.bed <- paste( c( out.dir, "lineage_genes-EUB.AFB.bed" ),
                          collapse = "" )

o.metadata.tsv <- paste( c( out.dir, "lineage_metadata-EUB.AFB.tsv" ),
                         collapse = "" )


if ( verb ) {
  cat( "Getting gene and TCC IDs and TCC expression file for \n" )
}

lineages <- unlist( strsplit( lineage, "," ))
gene.tcc.ids <- tibble::tibble( "trId" = character(), "geneId" = character() )


for ( lin in lineages ) {
  if ( verb ) {
    cat( c(lin, "lineage\n" ), sep = " " )
  }


  o.tcc.exp.tsv <- paste( c( out.dir, lin, "_tcc-exp-EUB.AFB.tsv.gz" ),
                        collapse = "" )

  i.tcc.exp <- paste( c( tcc.exp.dir, lin, tcc.exp.rds ), collapse = "" )
  i.gene.tcc.id <- paste( c( id.tcc.dir, lin, id.tcc.rds ), collapse = "" )


  current.ids <- get_gene_trans_ids_all( in_rds = i.gene.tcc.id )
  get_transcript_expression_file_all( in_rds = i.tcc.exp,
                                      ids = current.ids,
                                      out_tbl = o.tcc.exp.tsv )

  gene.tcc.ids <- gene.tcc.ids %>%
    dplyr::full_join( current.ids, by = colnames( gene.tcc.ids ))
}


if ( verb ) {
  cat( "Creating gene location file...\n" )
}

get_gene_location_file( in_gtf = annotation.gtf,
                        ids = gene.tcc.ids,
                        out_tbl = o.genes.loc.bed )

if ( verb ) {
  cat( "Creating PCs table...\n" )
}

pca.tbl <- get_pc_tbl( pca.path = pca.dir )


if ( verb ) {
  cat( "Creating metadata file...\n" )
}

get_metadata_file_all( in.tbl = metadata.tsv,
                       pca.tbl = pca.tbl,
                       out_tbl = o.metadata.tsv )

if ( verb ) {
  cat( "Done.\n" )
}
