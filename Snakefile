# PIPELINE ###################################################
# Decision:
# This pipeline will label sequences as Known, Potentially Novel, or Inconclusive.
# #############################################################

# INPUTS, METADATA & TRACEABILITY #############################

configfile: "config/config.yaml"

import csv

manifest = "data/manifest.tsv"

SAMPLES = []
READ1 = {}
READ2 = {}
SAMPLE_TYPE = {}

with open(manifest) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample = row["sample"]
        SAMPLES.append(sample)
        READ1[sample] = row["read1"]
        READ2[sample] = row["read2"]
        SAMPLE_TYPE[sample] = row.get("sample_type", "test")

rule all:
    input:
        expand("results/{sample}/call.json", sample=SAMPLES), # final decision artifact for a sample
        expand("results/{sample}/metrics.tsv", sample=SAMPLES) #monitoring / QA / dashboard-friendly table

## Step 1 : Raw read QC & Controls

rule S1_trim_filter:
    input:
        r1=lambda wc: READ1[wc.sample],
        r2=lambda wc: READ2[wc.sample]
    output:
        r1="results/{sample}/clean/{sample}_R1.fastq.gz",
        r2="results/{sample}/clean/{sample}_R2.fastq.gz"
    log:
        "logs/{sample}.cutadapt.log"
    conda: "workflow/envs/dada2.yaml"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.r1}) $(dirname {log})
        cutadapt -q 20 -m 100 \
          -o {output.r1} -p {output.r2} \
          {input.r1} {input.r2} > {log} 2>&1
        """ 

## Step 2 : Denoising / ASV inference + Taxonomic assignment
rule dada2:
    input:
        r1="results/{sample}/clean/{sample}_R1.fastq.gz",
        r2="results/{sample}/clean/{sample}_R2.fastq.gz",
        db="resources/db/silva_nr99_v138.2_toGenus_trainset.fa.gz",
        db_md5="resources/db/silva_nr99_v138.2_toGenus_trainset.fa.gz.md5"
    output:
        asv="results/{sample}/asv/asv_table.tsv",
        reps="results/{sample}/asv/rep_seqs.fasta",
        stats="results/{sample}/asv/denoise_stats.json",
        tax="results/{sample}/taxonomy/taxonomy.tsv"
    log:
        "logs/{sample}.dada2.log"
    threads : 1
    conda: "workflow/envs/dada2.yaml"
    params:
        # keep these simple; can be tuned later
        trunc_len_f=config.get("dada2", {}).get("trunc_len_f", 240),
        trunc_len_r=config.get("dada2", {}).get("trunc_len_r", 200)
    shell:
        r"""
        mkdir -p $(dirname {output.asv}) $(dirname {output.tax}) $(dirname {log})
        Rscript workflow/scripts/dada2_single_sample.R \
          --r1 {input.r1} --r2 {input.r2} \
          --db {input.db} --db_md5 {input.db_md5} \
          --trunc_len_f {params.trunc_len_f} --trunc_len_r {params.trunc_len_r} \
          --out_asv {output.asv} --out_reps {output.reps} --out_stats {output.stats} --out_tax {output.tax} \
          > {log} 2>&1
        """

rule db_checksum:
    input:
        "resources/db/silva_nr99_v138.2_toGenus_trainset.fa.gz"
    output:
        "resources/db/silva_nr99_v138.2_toGenus_trainset.fa.gz"
    shell:
        "md5sum {input} > {output}"

rule decide:
    input:
        tax="results/{sample}/taxonomy/taxonomy.tsv",
        asv="results/{sample}/asv/asv_table.tsv",
        stats="results/{sample}/asv/denoise_stats.json",
        manifest="data/manifest.tsv"
    output:
        call="results/{sample}/call.json",
        metrics="results/{sample}/metrics.tsv"
    log:
        "logs/{sample}.decide.log"
    conda: "workflow/envs/python.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.call}) $(dirname {log})
        python workflow/scripts/decide_novelty.py \
          --sample {wildcards.sample} \
          --manifest {input.manifest} \
          --tax {input.tax} --asv {input.asv} --stats {input.stats} \
          --out_call {output.call} --out_metrics {output.metrics} \
          > {log} 2>&1
        """







