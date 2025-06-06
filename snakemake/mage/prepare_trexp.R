#!/usr/bin/env Rscript

#### Prepare transcript expression file

## 1. Load libraries and arguments

library(optparse)
library(data.table)
# library(sQTLseekeR2)

option_list <- list(
    make_option(c("-g", "--group"), type = "character",
                help = "select sampleIds belonging to group", metavar = "CHARACTER"),
    make_option(c("-t", "--transcript_expr"), type = "character",
                help = "transcript expression file", metavar = "FILE"),
    make_option(c("-m", "--metadata"), type = "character",
                help = "file defining sample groups (sampleId, indId, group, [covariates])",
                metavar = "FILE"),
    make_option(c("-c", "--covariates"), action = "store_true",
                help = "prepare covariate file [default %default]",
                default = FALSE),
    make_option(c("-l", "--gene_location"), type = "character",
                help = "gene location file (chr, start, end, id)", metavar = "FILE"),
    make_option(c("-e", "--min_gene_expr"), type = "numeric", default = 1,
                help = "minimum gene expression [default %default]", metavar = "NUMERIC"),
    make_option(c("-p", "--min_proportion"), type = "numeric", default = 0.8,
                help = "minimum proportion of samples with gene expression above --min_gene_expr. [default %default]",
                metavar = "NUMERIC"),
    make_option(c("-i", "--min_transcript_expr"), type = "numeric", default = 0.1,
                help = "minimum transcript expression. [default %default]", metavar = "NUMERIC"),
    make_option(c("-d", "--min_dispersion"), type = "numeric", default = 0.1,
                help = "minimum dispersion of transcript relative expression. [default %default]",
                metavar = "NUMERIC"),
    make_option(c("--output_tre"), type = "character",
                help = "preprocessed transcript expression file", metavar = "FILE"),
    make_option(c("--output_gene"), type = "character",
                help = "updated gene location file", metavar = "FILE"),
    make_option(c("--output_cov"), type = "character",
                help = "prepared covariate file", metavar = "FILE"),
    make_option(c("-s", "--seed"), type = "numeric", help = "Set seed for random processess",
                metavar = "NUMERIC", default = 123),
    make_option(c("-v", "--verbose"), action = "store_true",
                help = "print genes and transcripts filtered out [default %default]",
                default = TRUE)
)

prepare.trans.exp <- function(te.df, min.transcript.exp = 0.01, min.gene.exp = 0.01,
                              min.prop = 0.8, min.dispersion = 0.1, verbose = FALSE)
{
    if(!all(c("geneId", "trId") %in% colnames(te.df))){
        stop("Missing column in 'te.df' : 'geneId' and 'trId' are required.")
    }
    if(min.transcript.exp < 0 || min.gene.exp < 0){
        stop("transcript/gene minimum expression should be greater or equal to 0.")
    }
    if(min.prop < 0 || min.prop > 1){
        stop("'min.prop' value should be in [0, 1].")
    }
    te.df$geneId <- as.character(te.df$geneId)
    te.df$trId <- as.character(te.df$trId)
    samples <- setdiff(colnames(te.df), c("chr", "start", "end", "geneId", "trId"))
    if(length(samples) < 5){
        stop("Not enough samples; at least 5 samples required (although at least 20 is recommended).")
    }
    if(length(samples) < 20){
        warning("Low sample size : it's recommended to have at least 20 samples.")
    }
    trans.to.keep <- apply(te.df[, samples], 1, function(r) any(r >= min.transcript.exp))
    if(all(!trans.to.keep)){
        stop("No transcripts with expression above threshold")
    }
    if(verbose && any(!trans.to.keep)){
        message("Filtered transcripts due to low expression : ",
                paste(te.df$trId[which(!trans.to.keep)], collapse = " "))
    }
    te.df <- te.df[which(trans.to.keep), ]
    nb.trans <- table(te.df$geneId)
    trans2 <- names(which(nb.trans > 1))
    if(verbose && any(nb.trans <= 1)){
        message("Filtered single-transcript genes : ",
                paste(setdiff(unique(te.df$geneId), trans2), collapse = " "))
    }
    te.df <- te.df[which(te.df$geneId %in% trans2), ]
    relativize.filter.dispersion <- function(df) {
        df[, samples] <- apply(df[, samples], 2, relativize, min.gene.exp = min.gene.exp)
        disp <- te.dispersion(df[, samples])
        non.na.p <- sum(!is.na(df[1, samples])) / length(samples)
        if (disp >= min.dispersion && non.na.p >= min.prop  &&
            nbDiffPt(df[, samples]) >= min(25, length(samples) * 0.8)){
            return(df)
        } else {
            return(data.frame())
        }
    }
    te.df <- plyr::ldply(lapply(unique(te.df$geneId), function(gene.i){
        df <- te.df[which(te.df$geneId == gene.i), ]
        relativize.filter.dispersion(df)
    }), identity)
    if(verbose && length(unique(te.df$geneId)) != length(trans2)){
        message("Filtered low exp/disp genes : ",
                paste(setdiff(trans2,unique(te.df$geneId)), collapse = " "))
    }
    if(nrow(te.df) == 0){
        stop("No genes found with suitable transcript expression.")
    }
    return(te.df)
}

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 2. Input files: transcript expression, sample groups, gene location

trans.expr.f <- opt$transcript_expr
metadata.f <- opt$metadata
genes.bed.f <- opt$gene_location
sel.group <- opt$group
out.tre.f <- opt$output_tre
out.gene.f <- opt$output_gene
out.cov.f <- opt$output_cov

if ( is.null(trans.expr.f) || is.null(metadata.f) || is.null(genes.bed.f) ||
     is.null(out.tre.f) || is.null (out.gene.f) ){
    print_help(opt_parser)
    stop("Missing/not found input files", call.= FALSE)
}

genes.bed <- read.table(genes.bed.f, header = TRUE, as.is = TRUE, sep = "\t")

## 3. Getting the IDs of the samples in the group of interest

metadata <- read.table(metadata.f, header = TRUE,
                       as.is = TRUE, sep = "\t")

subset.df <- subset(metadata, group == sel.group)                                   # Select samples of interest
subset.samples <- subset.df$sampleId

## 4. Prepare transcript expression

if( grepl("\\.gz$", trans.expr.f) ){
    trans.expr.f <- paste0("zcat < '", trans.expr.f, "'")
}

te.df <- as.data.frame(fread(input=trans.expr.f, header = TRUE, sep = "\t"))
colnames(te.df)[1:2] <- c("trId", "geneId")                                         # Proper 1,2 colnames
subset.samples <- subset.samples[subset.samples%in%colnames(te.df)]                 # Get samples that have quantifications
te.df <- te.df[, c("trId", "geneId", subset.samples)]                               # Select subgroup of samples = the ones from the group of interest
te.df <- subset(te.df, geneId %in% genes.bed$geneId)                                # Remove from te.df all genes that are not in genes.bed

set.seed(opt$seed)
tre.df <- prepare.trans.exp(te.df, min.gene.exp=opt$min_gene_expr,
                            min.transcript.exp=opt$min_transcript_expr,
                            min.dispersion=opt$min_dispersion,
                            min.prop=opt$min_proportion,
                            verbose=opt$verbose)                                  # Run

## 5. Prepare covariate file

if(opt$covariates) {
    if ( is.null(out.cov.f) ){
        print_help(opt_parser)
        stop("Missing output covariate file", call.= FALSE)
    }
    covariates.df <- subset.df[, setdiff(colnames(subset.df),
                                         c("sampleId", "indId", "group")), drop = FALSE]
    rownames(covariates.df) <- subset.df$indId
    all.na <- unlist(lapply(covariates.df, function(x){all(is.na(x))}))
    if (opt$verbose && sum(all.na) > 0){
        warning(sprintf("Covariates with NA values for all individuals were removed: %s",
                        paste(names(which(all.na)), collapse = ", ")))
    }
    covariates.df <- covariates.df[, !all.na, drop = FALSE]
    for(i in 1:ncol(covariates.df)){
        typ <- class(covariates.df[, i])
        if(typ == "character"){
            covariates.df[, i] <- as.factor(covariates.df[, i])
        } else if (typ == "numeric" || typ == "integer"){
            next
        } else {
            stop ("Covariates should be either 'numeric' or 'character'")
        }
    }
    types <- unlist(lapply(covariates.df, class))
    if (opt$verbose) {
        message("Covariate types:\n",
                paste(names(types), types, sep=": ", collapse = ", "))
    }
    save(covariates.df, file = out.cov.f)
} else {
    covariates.df <- NULL
    save(covariates.df, file = out.cov.f)
}

colnames(tre.df)[-c(1:2)] <- subset.df$indId                                        # Rename colnames to individual ID

## 6. Save result

save(tre.df, file = out.tre.f)                                                      # Save tre.df as RData

## 7. Preprocess for split

genes.bed <- subset(genes.bed, geneId %in% tre.df$geneId)                           # Remove from gene.bed all genes that are not in tre.df
write.table(genes.bed, file = out.gene.f, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t")                       # Write gene list

#### END
