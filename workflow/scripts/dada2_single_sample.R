#!/usr/bin/env Rscript

# Minimal DADA2 single-sample paired-end pipeline
# Inputs: trimmed/filtered FASTQ (paired-end) + SILVA toGenus trainset
# Outputs:
#   - asv_table.tsv
#   - rep_seqs.fasta
#   - denoise_stats.json
#   - taxonomy.tsv
#
# Notes:
# - Assumes cutadapt step already did primer/adaptor removal + basic quality filtering.
# - This script still performs DADA2 filtering/trimming (recommended) 

suppressPackageStartupMessages({
  library(dada2)
  library(jsonlite)
  library(optparse)
})
# -----------------------------
# STATS
# -----------------------------

stats <- list(
  sample_id=NULL,
  dada2_pipeline = "INCOMPLETE",
  fail_reason = NULL, 
  opt = NULL,

  step_filt_trim = NULL,
  reads_in=NULL,
  reads_after_filter=NULL,

  step_dada2 = NULL, 
  dada2_reads_input=NULL,
  dada2_reads_denoised = NULL,
  dada2_reads_merged= NULL,

  step_asv =NULL, 
  reads_nochim = NULL,
  num_asvs = NULL,

  params=list(
    trunc_len_f=NULL,
    trunc_len_r=NULL,
    max_ee_f=NULL,
    max_ee_r=NULL,
    trunc_q=NULL
  ),

  tools=list(
    dada2_version=as.character(packageVersion("dada2")),
    r_version=R.version.string
  ),

  db=list(
    name="SILVA", 
    version = "v138.2",
    trainset=NULL,
    md5=NULL

))

# -----------------------------
# Helpers
# -----------------------------

ensure_dir <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

read_md5 <- function(md5_path) {
  if (is.null(md5_path) || !file.exists(md5_path)) return(NA_character_)
  # md5sum output is "hash  filename"
  line <- readLines(md5_path, warn=FALSE)
  if (length(line) == 0) return(NA_character_)
  strsplit(line[1], "\\s+")[[1]][1]
}

# Try to infer sample_id from filename
infer_sample_id <- function(r1_path) {
  base <- basename(r1_path)
  # common patterns: SAMPLE_R1..., SAMPLE_1..., etc.
  base <- sub("\\.fastq\\.gz$", "", base)
  base <- sub("\\.fq\\.gz$", "", base)
  base <- sub("_R1.*$", "", base)
  base <- sub("_1.*$", "", base)
  base
}

graceful_exit <- function(
  stats
  ){
    # Ensure dir exsist.
    ensure_dir(stats$opt$out_asv); 
    ensure_dir(stats$opt$out_reps); 
    ensure_dir(stats$opt$out_stats); 
    ensure_dir(stats$opt$out_tax)

    # Empty ASV table
    empty_asv <- data.frame(
      asv_id = character(),
      sequence = character(),
      count = integer()
      )

    write.table(
      empty_asv,
      file = stats$opt$out_asv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    # Empty FASTA
    writeLines(character(), con=stats$opt$out_reps)

    # Empty taxonomy table
    tax_df <- data.frame(
      asv_id=character(),
      Kingdom=character(), 
      Phylum=character(), 
      Class=character(), 
      Order=character(),
      Family=character(), 
      Genus=character()
    )

    write.table(
      tax_df, 
      file=stats$opt$out_tax, 
      sep="\t", 
      quote=FALSE, 
      row.names=FALSE
      )

    write_json(stats, stats$opt$out_stats, pretty=TRUE, auto_unbox=TRUE)
    quit(save="no", status=0)
}

# -----------------------------
# CLI arguments
# -----------------------------
option_list <- list(
  make_option(c("--r1"), type="character", help="Path to R1 FASTQ.gz", metavar="FILE"),
  make_option(c("--r2"), type="character", help="Path to R2 FASTQ.gz", metavar="FILE"),
  make_option(c("--db"), type="character", help="Path to SILVA trainset .fa.gz", metavar="FILE"),
  make_option(c("--db_md5"), type="character", help="Path to db md5 file (optional)", metavar="FILE"),
  make_option(c("--trunc_len_f"), type="integer", default=240, help="Truncate length forward [default %default]"),
  make_option(c("--trunc_len_r"), type="integer", default=160, help="Truncate length reverse [default %default]"),
  make_option(c("--max_ee_f"), type="double", default=2, help="Max expected errors forward [default %default]"),
  make_option(c("--max_ee_r"), type="double", default=2, help="Max expected errors reverse [default %default]"),
  make_option(c("--trunc_q"), type="integer", default=2, help="Truncate at quality <= truncQ [default %default]"),
  make_option(c("--out_asv"), type="character", help="Output ASV table TSV", metavar="FILE"),
  make_option(c("--out_reps"), type="character", help="Output representative sequences FASTA", metavar="FILE"),
  make_option(c("--out_stats"), type="character", help="Output denoise stats JSON", metavar="FILE"),
  make_option(c("--out_tax"), type="character", help="Output taxonomy TSV", metavar="FILE")
)

opt <- parse_args(OptionParser(option_list = option_list))

`%||%` <- function(a, b) if (!is.null(a)) a else b

required <- c(
  "r1", "r2", "db",
  "out_asv", "out_reps", "out_stats", "out_tax"
)

missing <- required[
  !nzchar(sapply(required, function(x) opt[[x]] %||% ""))
]

if (length(missing) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(missing, collapse = ", ")
  )
}

stats$opt <- opt

stats$params$trunc_len_f=opt$trunc_len_f
stats$params$trunc_len_r=opt$trunc_len_r
stats$params$max_ee_f=opt$max_ee_f
stats$params$max_ee_r=opt$max_ee_r
stats$params$trunc_q=opt$trunc_q

stats$db$trainset=opt$db
stats$db$md5=read_md5(opt$db_md5)

sample_id <- infer_sample_id(opt$r1)
stats$sample_id <- infer_sample_id(opt$r1)

# -----------------------------
# DADA2 pipeline
# -----------------------------
# 1) Filter and trim into temp files (keeps this script self-contained)
tmpdir <- tempfile(paste0("dada2_", sample_id, "_"))
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)

filtF <- file.path(tmpdir, paste0(sample_id, "_filt_F.fastq.gz"))
filtR <- file.path(tmpdir, paste0(sample_id, "_filt_R.fastq.gz"))

# Filtering - filters and trims based on user defined criteria
# Output is a matrix, row = processed sample or file pair, column = reads.in(raw), read.out(processes) & compressed filterd reads fastq. 
filt_out <- filterAndTrim(
  fwd=opt$r1, filt=filtF,
  rev=opt$r2, filt.rev=filtR,
  truncLen=c(opt$trunc_len_f, opt$trunc_len_r), # reads shorter than this are discarted
  maxEE=c(opt$max_ee_f, opt$max_ee_r), # MaximumExpectedErrors, read with more than the threshold of expected errors will be discarted
  truncQ=opt$trunc_q, # Truncate reads at the 1st instance the quality is <= Q
  rm.phix=TRUE, # Removes reads matching to PhiX Bacteriophage
  compress=TRUE, # Compresses output
  multithread=FALSE 
)

reads_in <- as.integer(filt_out[1, "reads.in"])
reads_out <- as.integer(filt_out[1, "reads.out"])

stats$reads_in <- reads_in
stats$reads_after_filter <- reads_out

# If no reads survive, write minimal outputs and exit gracefully
if (is.na(reads_out) || reads_out == 0) {
  stats$step_filt_trim <- "FAIL"
  stats$fail_reason <- "No reads after filtering; outputs are empty."
  graceful_exit(stats=stats)
} else {
  stats$step_filt_trim <- "PASS"
}

# 2) Learn errors
# model learns errors by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution
errF <- learnErrors(filtF, multithread=TRUE)
errR <- learnErrors(filtR, multithread=TRUE)

# 3) Dereplicate
# dereplicating amplicon sequences
derepF <- derepFastq(filtF)
derepR <- derepFastq(filtR)

# 4) DADA 
# inference to correct for Illumina sequenced amplicon errors. Actual Denoising.
dadaF <- dada(derepF, err=errF, multithread=TRUE)
dadaR <- dada(derepR, err=errR, multithread=TRUE)

# 5) Merge pairs
# merge each denoised pair of forward and reverse reads, rejecting any pairs which do not sufficiently overlap or which contain too many mismatches
mergers <- mergePairs(
    dadaF, derepF, dadaR, derepR, verbose = FALSE)

# Track counts (single sample)
getN <- function(x) sum(getUniques(x))

stats$dada2_reads_input <- reads_out
stats$dada2_reads_denoised <- getN(dadaF)
stats$dada2_reads_merged <- sum(mergers$abundance)
stats$num_merged_variants <- nrow(mergers)


if (nrow(mergers) == 0L) {
  stats$step_dada2 <- "FAIL"
  stats$fail_reason <-"No forward and reverse reads could be merged."
  graceful_exit ( stats=stats  )
} else {
  stats$step_dada2 <- "PASS"
}


# 6) Make sequence table + remove chimeras
# Creates ASV table
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

stats$reads_nochim <- sum(seqtab.nochim)
stats$num_asvs <- ncol(seqtab.nochim)

if (stats$num_asvs == 0L) {
  stats$step_asv <- "FAIL"
  stats$fail_reason <- "No asvs after chimera removal."
  graceful_exit(stats)
} else {
  stats$step_asv <- "PASS"
}


# 7) Taxonomy assignment (to genus using SILVA trainset)
tax <- assignTaxonomy(seqtab.nochim, opt$db, multithread=FALSE)

# -----------------------------
# Write outputs
# -----------------------------
# Makesure the dada2 pipeline files are exsist
ensure_dir(opt$out_asv)
ensure_dir(opt$out_reps)
ensure_dir(opt$out_stats)
ensure_dir(opt$out_tax)

# 1) Make ASV Table
# ASV IDs: stable names derived from sequence hashes (short)
# Sets up column (seq) and rows (asv_ids) of ASV table
seqs <- colnames(seqtab.nochim)
asv_ids <- paste0("ASV", seq_along(seqs))

# ASV table: one sample row (wide) OR long (I recommend long for simplicity)
# We'll write long: asv_id, sequence, count
counts <- as.integer(seqtab.nochim[1, ])
asv_table <- data.frame(
  asv_id=asv_ids,
  sequence=seqs,
  count=counts,
  stringsAsFactors=FALSE
)
write.table(asv_table, file=opt$out_asv, sep="\t", quote=FALSE, row.names=FALSE)

# 2) Make representative sequences FASTA
fasta_lines <- c(rbind(paste0(">", asv_ids), seqs))
writeLines(fasta_lines, con=opt$out_reps)

# 3) Taxonomy TSV
# tax is a matrix with columns (Kingdom..Genus depending)
tax_mat <- as.matrix(tax)
# Ensure columns exist
wanted_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus")
for (cname in wanted_cols) {
  if (!(cname %in% colnames(tax_mat))) {
    tax_mat <- cbind(tax_mat, setNames(matrix(NA_character_, nrow=nrow(tax_mat), ncol=1), cname))
  }
}
# make dataframe
tax_out <- data.frame(
  asv_id=asv_ids,
  tax_mat[, wanted_cols, drop=FALSE],
  stringsAsFactors=FALSE
)
write.table(tax_out, file=opt$out_tax, sep="\t", quote=FALSE, row.names=FALSE)

stats$dada2_pipeline <- "COMPLETE"
write_json(stats, opt$out_stats, pretty=TRUE, auto_unbox=TRUE)

message(
  "Done. Sample: ",
  stats$sample_id,
  " | ASVs: ",
  stats$num_asvs,
  " | reads_nochim: ",
  stats$reads_nochim
)
