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

opt <- parse_args(OptionParser(option_list=option_list))

`%||%` <- function(a, b) if (!is.null(a)) a else b

required <- c("r1","r2","db","out_asv","out_reps","out_stats","out_tax")
missing <- required[!nzchar(sapply(required, function(x) opt[[x]] %||% ""))]
if (length(missing) > 0) {
  stop(paste0("Missing required args: ", paste(missing, collapse=", ")))
}



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

sample_id <- infer_sample_id(opt$r1)

# -----------------------------
# DADA2 pipeline
# -----------------------------
# 1) Filter and trim into temp files (keeps this script self-contained)
tmpdir <- tempfile(paste0("dada2_", sample_id, "_"))
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)

filtF <- file.path(tmpdir, paste0(sample_id, "_filt_F.fastq.gz"))
filtR <- file.path(tmpdir, paste0(sample_id, "_filt_R.fastq.gz"))

# Filtering
filt_out <- filterAndTrim(
  fwd=opt$r1, filt=filtF,
  rev=opt$r2, filt.rev=filtR,
  truncLen=c(opt$trunc_len_f, opt$trunc_len_r),
  maxEE=c(opt$max_ee_f, opt$max_ee_r),
  truncQ=opt$trunc_q,
  rm.phix=TRUE,
  compress=TRUE,
  multithread=FALSE
)

reads_in <- as.integer(filt_out[1, "reads.in"])
reads_out <- as.integer(filt_out[1, "reads.out"])

# If no reads survive, write minimal outputs and exit gracefully
if (length(missing) > 0) {
  stop(paste0("Missing required args: ", paste(missing, collapse=", ")))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

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

sample_id <- infer_sample_id(opt$r1)

# -----------------------------
# DADA2 pipeline
# -----------------------------
# 1) Filter and trim into temp files (keeps this script self-contained)
tmpdir <- tempfile(paste0("dada2_", sample_id, "_"))
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)

filtF <- file.path(tmpdir, paste0(sample_id, "_filt_F.fastq.gz"))
filtR <- file.path(tmpdir, paste0(sample_id, "_filt_R.fastq.gz"))

# Filtering
filt_out <- filterAndTrim(
  fwd=opt$r1, filt=filtF,
  rev=opt$r2, filt.rev=filtR,
  truncLen=c(opt$trunc_len_f, opt$trunc_len_r),
  maxEE=c(opt$max_ee_f, opt$max_ee_r),
  truncQ=opt$trunc_q,
  rm.phix=TRUE,
  compress=TRUE,
  multithread=TRUE
)

reads_in <- as.integer(filt_out[1, "reads.in"])
reads_out <- as.integer(filt_out[1, "reads.out"])

# If no reads survive, write minimal outputs and exit gracefully
if (is.na(reads_out) || reads_out == 0) {
  ensure_dir(opt$out_asv); ensure_dir(opt$out_reps); ensure_dir(opt$out_stats); ensure_dir(opt$out_tax)

  # Empty ASV table
  write.table(data.frame(), file=opt$out_asv, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  # Empty FASTA
  writeLines(character(), con=opt$out_reps)

  # Empty taxonomy table
  tax_df <- data.frame(
    asv_id=character(),
    Kingdom=character(), Phylum=character(), Class=character(), Order=character(),
    Family=character(), Genus=character()
  )
  write.table(tax_df, file=opt$out_tax, sep="\t", quote=FALSE, row.names=FALSE)

  stats <- list(
    sample_id=sample_id,
    reads_in=reads_in,
    reads_after_filter=reads_out,
    dada2_reads_input=reads_out,
    dada2_reads_merged=0,
    reads_nochim=0,
    num_asvs=0,
    db=list(name="SILVA", trainset=opt$db, md5=read_md5(opt$db_md5)),
    note="No reads after filtering; outputs are empty."
  )
  write_json(stats, opt$out_stats, pretty=TRUE, auto_unbox=TRUE)
  quit(save="no", status=0)
}

# 2) Learn errors
errF <- learnErrors(filtF, multithread=TRUE)
errR <- learnErrors(filtR, multithread=TRUE)

# 3) Dereplicate
derepF <- derepFastq(filtF)
derepR <- derepFastq(filtR)
names(derepF) <- sample_id
names(derepR) <- sample_id

# 4) DADA inference
dadaF <- dada(derepF, err=errF, multithread=TRUE)
dadaR <- dada(derepR, err=errR, multithread=TRUE)

# 5) Merge pairs
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose=FALSE)

# 6) Make sequence table + remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

# Track counts (single sample)
getN <- function(x) sum(getUniques(x))
dada2_reads_input <- reads_out
dada2_reads_denoised <- getN(dadaF) # forward uniques count proxy
dada2_reads_merged <- sum(seqtab)
reads_nochim <- sum(seqtab.nochim)
num_asvs <- ncol(seqtab.nochim)

# 7) Taxonomy assignment (to genus using SILVA trainset)
tax <- assignTaxonomy(seqtab.nochim, opt$db, multithread=FALSE)

# -----------------------------
# Write outputs
# -----------------------------
ensure_dir(opt$out_asv)
ensure_dir(opt$out_reps)
ensure_dir(opt$out_stats)
ensure_dir(opt$out_tax)

# ASV IDs: stable names derived from sequence hashes (short)
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

# Representative sequences FASTA
fasta_lines <- c(rbind(paste0(">", asv_ids), seqs))
writeLines(fasta_lines, con=opt$out_reps)

# Taxonomy TSV
# tax is a matrix with columns (Kingdom..Genus depending)
tax_mat <- as.matrix(tax)
# Ensure columns exist
wanted_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus")
for (cname in wanted_cols) {
  if (!(cname %in% colnames(tax_mat))) {
    tax_mat <- cbind(tax_mat, setNames(matrix(NA_character_, nrow=nrow(tax_mat), ncol=1), cname))
  }
}
tax_out <- data.frame(
  asv_id=asv_ids,
  tax_mat[, wanted_cols, drop=FALSE],
  stringsAsFactors=FALSE
)
write.table(tax_out, file=opt$out_tax, sep="\t", quote=FALSE, row.names=FALSE)

# Stats JSON (provenance-ish)
stats <- list(
  sample_id=sample_id,
  reads_in=reads_in,
  reads_after_filter=reads_out,
  dada2_reads_input=dada2_reads_input,
  dada2_reads_merged=dada2_reads_merged,
  reads_nochim=reads_nochim,
  num_asvs=num_asvs,
  params=list(
    trunc_len_f=opt$trunc_len_f,
    trunc_len_r=opt$trunc_len_r,
    max_ee_f=opt$max_ee_f,
    max_ee_r=opt$max_ee_r,
    trunc_q=opt$trunc_q
  ),
  tools=list(
    dada2_version=as.character(packageVersion("dada2")),
    r_version=R.version.string
  ),
  db=list(
    name="SILVA",
    version="v138.2_toGenus_trainset",
    trainset=opt$db,
    md5=read_md5(opt$db_md5)
  )
)
write_json(stats, opt$out_stats, pretty=TRUE, auto_unbox=TRUE)

message("Done. Sample: ", sample_id, " | ASVs: ", num_asvs, " | reads_nochim: ", reads_nochim)
