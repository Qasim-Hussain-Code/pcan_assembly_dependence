"""Contiguity statistics for every cached assembly in a manifest.

Contig N50 here splits sequences at every run of N, the same rule used for the
fragmented replicates of arm 1, so that published and perturbed assemblies sit
on one axis. NCBI's contig N50 is copied beside it from the dataset report.
"""
import argparse
import json
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import padlib as P  # noqa: E402


def sequence_roles(path):
    roles = {}
    if not os.path.exists(path):
        return roles
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            for key in ("genbankAccession", "refseqAccession"):
                if r.get(key):
                    roles[r[key]] = r
    return roles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    man = pd.read_csv(args.manifest, sep="\t", dtype=str)
    rows = []
    for acc in man["accession"]:
        base = P.data_dir("assemblies", acc)
        recs = P.read_fasta(base + ".fna.gz")
        roles = sequence_roles(base + ".sequence_report.jsonl")
        with open(base + ".report.json") as fh:
            rep = [r for r in json.load(fh)["reports"] if r["accession"] == acc][0]
        seq_lengths, contig_lengths = [], []
        gap_runs = gap_bases = at = acgt = 0
        nuclear = mito = 0
        for name, seq in recs:
            seq_lengths.append(len(seq))
            contig_lengths.extend(P.contig_lengths_split_at_gaps(seq))
            gap_bases += seq.upper().count("N")
            gap_runs += len(re.findall("[Nn]+", seq))
            s = seq.upper()
            at += s.count("A") + s.count("T")
            acgt += sum(s.count(b) for b in "ACGT")
            r = roles.get(name, {})
            if r.get("role") == "assembled-molecule":
                if r.get("assignedMoleculeLocationType") == "Chromosome":
                    nuclear += 1
                elif r.get("assignedMoleculeLocationType") == "Mitochondrion":
                    mito += 1
        stats = rep.get("assembly_stats", {})
        rows.append({
            "accession": acc,
            "organism": rep.get("organism", {}).get("organism_name", ""),
            "assembly_level": rep.get("assembly_info", {}).get("assembly_level", ""),
            "n_sequences": len(seq_lengths),
            "total_length": sum(seq_lengths),
            "sequence_n50": P.nx(seq_lengths),
            "contig_n50": P.nx(contig_lengths),
            "contig_l50": P.lx(contig_lengths),
            "n_contigs": len(contig_lengths),
            "n_gap_runs": gap_runs,
            "gap_bases": gap_bases,
            "nuclear_chromosome_sequences": nuclear,
            "mitochondrial_sequences": mito,
            "ncbi_contig_n50": stats.get("contig_n50", ""),
            "ncbi_scaffold_n50": stats.get("scaffold_n50", ""),
            "genome_at": round(at / acgt, 4) if acgt else "",
        })
        P.log("%s: %d sequences, contig N50 %d" % (acc, len(seq_lengths), rows[-1]["contig_n50"]))
    P.atomic_write_tsv(pd.DataFrame(rows), args.out)


if __name__ == "__main__":
    main()
