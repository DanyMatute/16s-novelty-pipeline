#!/usr/bin/env python3
import argparse
import csv
import json
import os
from datetime import datetime

def ensure_dir(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)

def read_manifest(manifest_path: str) -> dict:
    """
    Returns dict: sample -> row dict
    Expects columns: sample, read1, read2, optional sample_type
    """
    m = {}
    with open(manifest_path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            m[row["sample"]] = row
    return m

def load_json(path: str) -> dict:
    with open(path) as f:
        return json.load(f)

def read_asv_table(asv_path: str) -> list[dict]:
    """
    Expects DADA2 script output (long format):
    asv_id, sequence, count
    """
    rows = []
    with open(asv_path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            # tolerate empty tables
            if not row:
                continue
            # counts might be empty for empty file
            if "count" in row and row["count"] != "":
                row["count"] = int(float(row["count"]))
            rows.append(row)
    return rows

def read_taxonomy(tax_path: str) -> list[dict]:
    """
    Expects taxonomy.tsv with:
    asv_id, Kingdom, Phylum, Class, Order, Family, Genus
    """
    rows = []
    with open(tax_path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if not row:
                continue
            rows.append(row)
    return rows

def main():

    # 1) Set Up Arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--tax", required=True)
    ap.add_argument("--asv", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--out_call", required=True)
    ap.add_argument("--out_metrics", required=True)

    # decision thresholds (simple defaults; tune later)
    ap.add_argument("--min_reads_nochim", type=int, default=200)
    ap.add_argument("--min_asvs", type=int, default=5)
    ap.add_argument("--novel_genus_label", default="NA")  # placeholder for future
    args = ap.parse_args()

    # 2) Parce Arguments & Get sample type
    # Manifest
    manifest = read_manifest(args.manifest)
    if args.sample not in manifest:
        raise SystemExit(f"Sample '{args.sample}' not found in manifest: {args.manifest}")

    # Sample type
    sample_type = (manifest[args.sample].get("sample_type") or "test").strip().lower()

    # Denoising Stats
    stats = load_json(args.stats) if os.path.exists(args.stats) else {}
    reads_in = int(stats.get("reads_in", 0) or 0)
    reads_after_filter = int(stats.get("reads_after_filter", 0) or 0)
    reads_nochim = int(stats.get("reads_nochim", 0) or 0)
    num_asvs = int(stats.get("num_asvs", 0) or 0)

    # ASV & Taxonomy
    asvs = read_asv_table(args.asv) if os.path.exists(args.asv) else []
    tax = read_taxonomy(args.tax) if os.path.exists(args.tax) else []

    # 3) Simple “evidence” summaries
    # Most Aboundant ASVs
    top_asvs = sorted(asvs, key=lambda r: r.get("count", 0), reverse=True)[:5]
    # count genus assignments (ignore blanks)
    genus_calls = [r.get("Genus", "") for r in tax if (r.get("Genus") or "").strip() not in ("", "NA", "NaN")]
    num_genus_assigned = len(genus_calls)

    # 4) ---- Decision logic (workable v0) ----
    # Rule 0: If too few reads, can't conclude anything
    reasons = [] # collects reasons supporting decisions.
    call = "INCONCLUSIVE" # Automatically Inconclusive

    if sample_type == "negative":
        # For a negative control, any substantial signal is suspicious as it should be a blank.
        # if the number of reads after chimera removal & number of ASV is bigger than the threshold then the negative is not blank
        if reads_nochim >= args.min_reads_nochim and num_asvs >= args.min_asvs:
            call = "INCONCLUSIVE"
            reasons.append("negative_control_has_signal_investigate_contamination")
        # if the number of reads after chimera removal & number of ASV is smaller than the threshold then the negative is mostly blank, is therefore inconclusive as its novelty cant be classified, but this is expected and hoped for. 
        else:
            call = "INCONCLUSIVE"
            reasons.append("negative_control_low_signal_expected")
    else:
        # Non-negative samples
        # If the number of reads after chimeral removal is less than the threshold then many reads were amplicon errors caught by dada2.
        if reads_nochim < args.min_reads_nochim:
            call = "INCONCLUSIVE"
            reasons.append("low_reads_after_dada2_nochim")
        # if number of asc is less than threshold, then we can not be confident that the tanaomic identity was correctly derermined. We need many asv hits to a taxanomic identity. 
        elif num_asvs < args.min_asvs:
            call = "INCONCLUSIVE"
            reasons.append("too_few_asvs_for_confident_taxonomy")
        elif num_genus_assigned == 0:
            call = "POTENTIALLY_NOVEL"
            reasons.append("no_genus_assignment_from_silva_trainset")\
        # If there are enough reads_nochim, num_asvs, and a genus was assigned then the sample is known in the db. 
        else:
            call = "KNOWN"
            reasons.append("genus_assigned_using_silva_trainset")

    # Candidate novelty list (very conservative in v0)
    novel_candidates = []
    if call == "POTENTIALLY_NOVEL":
        for r in top_asvs:
            novel_candidates.append({
                "asv_id": r.get("asv_id"),
                "count": r.get("count", 0),
                "note": "no_genus_assignment_v0"
            })

    call_obj = {
        "sample_id": args.sample,
        "sample_type": sample_type,
        "call": call,
        "reasons": reasons,
        "novel_candidates": novel_candidates,
        "summary": {
            "reads_in": reads_in,
            "reads_after_filter": reads_after_filter,
            "reads_nochim": reads_nochim,
            "num_asvs": num_asvs,
            "num_genus_assigned": num_genus_assigned
        },
        "provenance": {
            "created_at": datetime.utcnow().isoformat() + "Z",
            "pipeline_version": "v0-demo",
            "note": "screening call; novelty indicates candidates for follow-up, not formal species description"
        }
    }

    # Metrics TSV (one row)
    metrics_header = [
        "sample_id","sample_type","call",
        "reads_in","reads_after_filter","reads_nochim",
        "num_asvs","num_genus_assigned",
        "top_asv_ids"
    ]
    top_ids = ",".join([r.get("asv_id","") for r in top_asvs if r.get("asv_id")])

    ensure_dir(args.out_call)
    with open(args.out_call, "w") as f:
        json.dump(call_obj, f, indent=2)

    ensure_dir(args.out_metrics)
    with open(args.out_metrics, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(metrics_header)
        w.writerow([
            args.sample, sample_type, call,
            reads_in, reads_after_filter, reads_nochim,
            num_asvs, num_genus_assigned,
            top_ids
        ])

if __name__ == "__main__":
    main()
