"""Arms 2 and 3: build the strain pairs (one long-read and one short-read
assembly of the same S. cerevisiae strain) and the Illumina runs for arm 3,
from the primary sources, by the rules in config/analysis_plan.md.

Sources, all fetched here and logged in logs/downloads.tsv:
  O'Donnell et al. 2023 (ScRAP), Supplementary Table 1: long-read assemblies
    and their standardized names
  Peter et al. 2018, Supplementary Tables S1 and S17: standardized names,
    isolate names, ploidy, zygosity, and the reads file behind each assembly
  the 1002 Yeast Genomes archive of Peter et al. 2018 assemblies (one 4.0 GB
    tar.gz, checked against the project's md5.txt)
  NCBI Datasets: status, technology and contig N50 of every ScRAP accession
  ENA: the file report of study PRJEB13017, Peter et al. 2018's reads

Usage:
    python scripts/09_build_pairs.py [--force]

Writes config/arm2_pairs.tsv, config/arm3_runs.tsv, results/arm2/pair_candidates.tsv,
and the extracted short-read assemblies in data/assemblies/peter2018/.
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import padlib as P  # noqa: E402

SCRAP_XLSX = ("scrap2023_supplementary_tables.xlsx",
              "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01459-y/MediaObjects/41588_2023_1459_MOESM4_ESM.xlsx")
PETER_XLS = ("peter2018_supplementary_tables.xls",
             "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0030-5/MediaObjects/41586_2018_30_MOESM3_ESM.xls")
PETER_ARCHIVE = ("1011Assemblies.tar.gz", "http://1002genomes.u-strasbg.fr/files/1011Assemblies.tar.gz")
PETER_MD5 = "http://1002genomes.u-strasbg.fr/files/md5.txt"
ENA_STUDY = "PRJEB13017"
ENA_FIELDS = ("run_accession,sample_accession,sample_title,library_name,instrument_model,library_layout,"
              "read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp")
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/%s/dataset_report?filters.assembly_version=all_assemblies"


def tables(*parts):
    return P.data_dir("tables", *parts)


def log_download(path, url, expect=""):
    sha = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)
    ok = "not published"
    if expect:
        md5 = hashlib.md5()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                md5.update(chunk)
        ok = "yes" if md5.hexdigest() == expect else "no"
    row = [P.rel(path), url, str(os.path.getsize(path)), sha.hexdigest(), ("md5:" + expect) if expect else "none",
           ok, "not given", P.utc_now()]
    with open(P.repo("logs", "downloads.tsv"), "a") as fh:
        fh.write("\t".join(row) + "\n")
    return ok


def fetch(name, url, sub="", expect=""):
    dest = tables(sub, name) if sub else tables(name)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if not os.path.exists(dest):
        P.log("downloading %s" % url)
        subprocess.run(["curl", "-fsSL", "--retry", "5", "--retry-delay", "15", "-o", dest + ".part", url], check=True)
        os.replace(dest + ".part", dest)
        if log_download(dest, url, expect) == "no":
            os.remove(dest)
            sys.exit("%s: md5 does not match %s" % (name, expect))
    return dest


def get_json(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:      # network errors are retried, then raised
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def norm_name(s):
    return re.sub(r"[\s_\-.]", "", str(s).lower())


def scrap_table(path):
    t = pd.read_excel(path, sheet_name="Supp Table 1", header=2)
    t = t[t["standardized_name"].notna()].copy()
    t["accession"] = t["assembly_accession"].astype(str).str.extract(r"(GC[AF]_\d+\.\d+)")[0]
    return t


def peter_tables(path):
    s1 = pd.read_excel(path, sheet_name="Table S1", header=3)
    s1 = s1[s1["Standardized name"].notna()].copy()
    s1["code"] = s1["Standardized name"].astype(str).str.strip()
    s17 = pd.read_excel(path, sheet_name="Table S17", header=2).iloc[:, :5]
    s17.columns = ["isolate", "reads_file", "n50", "total_bases", "total_contigs"]
    s17 = s17[s17["isolate"].notna()]
    return s1, s17


def ncbi_reports(accs):
    """Dataset report per accession, cached as JSON lines."""
    cache = tables("scrap_ncbi_reports.jsonl")
    have = {}
    if os.path.exists(cache):
        with open(cache) as fh:
            for line in fh:
                r = json.loads(line)
                have[r["queried"]] = r
    with open(cache, "a") as fh:
        for a in accs:
            if a in have:
                continue
            d = get_json(NCBI % a)
            rep = [r for r in d.get("reports", []) if r.get("accession") == a]
            rec = {"queried": a, "report": rep[0] if rep else None, "messages": d.get("messages")}
            fh.write(json.dumps(rec) + "\n")
            have[a] = rec
            time.sleep(0.4)
    return have


def archive_members(archive):
    listing = tables("peter2018", "members.txt")
    if not os.path.exists(listing):
        P.log("listing the Peter et al. 2018 archive")
        out = subprocess.run(["tar", "-tzf", archive], check=True, stdout=subprocess.PIPE,
                             universal_newlines=True).stdout
        P.atomic_write_text(listing, out)
    with open(listing) as fh:
        return [l.strip() for l in fh if l.strip()]


def ena_runs():
    path = tables("ena_%s_runs.tsv" % ENA_STUDY)
    if not os.path.exists(path):
        url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=%s&result=read_run&fields=%s&format=tsv"
               % (ENA_STUDY, ENA_FIELDS))
        subprocess.run(["curl", "-fsSL", "--retry", "5", "-o", path + ".part", url], check=True)
        os.replace(path + ".part", path)
        log_download(path, url)
    return pd.read_csv(path, sep="\t", dtype=str)


def submitted_stem(name):
    """BCM_AABOSW_6_1_C2420ACXX.IND2_clean.fastq.gz -> BCM_AABOSW_6_C2420ACXX.IND2,
    the form Table S17 writes: the mate number is dropped."""
    base = os.path.basename(name).replace("_clean.fastq.gz", "").replace(".fastq.gz", "")
    return re.sub(r"^(BCM_[A-Za-z0-9]+?OSW_\d+)_[12]_", r"\1_", base)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if P.is_done("build_pairs") and not args.force:
        P.log("pairs already built, skipping (use --force to redo)")
        return

    scrap = scrap_table(fetch(*SCRAP_XLSX))
    s1, s17 = peter_tables(fetch(*PETER_XLS))
    # md5.txt is in the BSD layout: MD5 (<file>) = <hash>
    md5_lines = urllib.request.urlopen(PETER_MD5, timeout=60).read().decode().splitlines()
    expect = [m.group(2) for m in (re.match(r"MD5 \((.+)\) = ([0-9a-f]{32})", l.strip()) for l in md5_lines)
              if m and m.group(1) == PETER_ARCHIVE[0]]
    if not expect:
        sys.exit("no md5 for %s on the server" % PETER_ARCHIVE[0])
    archive = fetch(*PETER_ARCHIVE, sub="peter2018", expect=expect[0])
    members = archive_members(archive)
    sd6 = pd.read_csv(tables("supp_data6_assemblies.tsv"), sep="\t", dtype=str)
    sd6_peter = set(sd6.loc[sd6["Dataset"] == "Peter et al. genomes", "Assembly"].str.strip())
    sd6_ncbi = set(sd6.loc[sd6["Dataset"] == "NCBI assemblies", "Assembly"].str.strip())

    peter = s1.set_index("code")
    reports = ncbi_reports(sorted(scrap["accession"].dropna().unique()))
    cand_rows, pairs = [], []
    for (code, strain), g in scrap.groupby(["standardized_name", "strain_name"], sort=True):
        code = str(code).strip()
        row = {"code": code, "scrap_strain_name": strain}
        if code not in peter.index:
            row.update(decision="dropped", reason="standardized name not in Peter et al. 2018 Table S1")
            cand_rows.append(row)
            continue
        p = peter.loc[code]
        row["peter_isolate_name"] = p["Isolate name"]
        if norm_name(strain) != norm_name(p["Isolate name"]):
            row.update(decision="dropped",
                       reason="ScRAP strain name differs from the Peter et al. 2018 isolate name (a derivative "
                              "or another strain under the same code)")
            cand_rows.append(row)
            continue
        options = []
        for _, a in g.iterrows():
            acc = a["accession"]
            rep = reports.get(acc, {}).get("report") if isinstance(acc, str) else None
            if not rep:
                options.append((acc, a["assembly"], a["assembly.type"], None, "not served by NCBI Datasets"))
                continue
            info, st = rep.get("assembly_info", {}), rep.get("assembly_stats", {})
            tech = info.get("sequencing_tech", "")
            if info.get("assembly_status") != "current":
                why = "NCBI status %s" % info.get("assembly_status")
            elif not re.search(r"nanopore|pacbio|\bONT\b", tech, re.I):
                why = "technology not long-read: %s" % tech
            else:
                why = ""
            options.append((acc, a["assembly"], a["assembly.type"], int(st.get("contig_n50", 0) or 0), why, tech,
                            int(st.get("total_sequence_length", 0) or 0)))
        usable = [o for o in options if len(o) > 5 and not o[4]]
        row["long_read_options"] = "; ".join("%s %s (%s)%s" % (o[0], o[1], o[2], (" contig N50 %d" % o[3]) if o[3] else "")
                                             + ((" [" + o[4] + "]") if o[4] else "") for o in options)
        if not usable:
            row.update(decision="dropped", reason="no long-read assembly that NCBI serves as current")
            cand_rows.append(row)
            continue
        best = sorted(usable, key=lambda o: (-o[3], o[0]))[0]
        # Members are named <code>_<lane>.re.fa for most strains, <code>.re.fa for
        # the externally sequenced ones, and in one case <code>_.re.fa; all
        # three name the code unambiguously.
        stem = re.escape(code.replace("SACE_", ""))
        hits = [m for m in members if re.match(r"^GENOMES_ASSEMBLED/%s(_\d*)?\.re\.fa$" % stem, m)]
        if len(hits) != 1:
            row.update(decision="dropped", reason="%d members of the Peter et al. 2018 archive match the code" % len(hits))
            cand_rows.append(row)
            continue
        row.update(decision="paired", reason="")
        cand_rows.append(row)
        pairs.append({
            "code": code, "isolate_name": p["Isolate name"], "long_read_accession": best[0],
            "long_read_assembly_name": best[1], "long_read_assembly_type": best[2], "long_read_tech": best[5],
            "long_read_contig_n50_ncbi": best[3], "long_read_length": best[6],
            "long_read_in_sd6": "yes" if best[0] in sd6_ncbi else "no",
            "short_read_source": "Peter et al. 2018, 1011Assemblies.tar.gz", "short_read_member": hits[0],
            "short_read_in_sd6": "yes" if code.replace("SACE_", "") in sd6_peter else "no",
            "peter_ploidy": p["Ploidy"], "peter_zygosity": p["Zygosity"], "peter_aneuploidies": p["Aneuploidies"],
            "scrap_sequenced_as": g.iloc[0]["sequenced_as"], "scrap_zygosity": g.iloc[0]["original_Zygosity"],
            "scrap_ploidy": g.iloc[0]["ploidy"],
            "match_evidence": "standardized name %s in ScRAP Supplementary Table 1 and Peter et al. 2018 Table S1; "
                              "isolate name %s in both" % (code, p["Isolate name"])})
    cand = pd.DataFrame(cand_rows)
    pairs = pd.DataFrame(pairs).sort_values("code")
    P.atomic_write_tsv(cand, P.repo("results", "arm2", "pair_candidates.tsv"))
    P.atomic_write_tsv(pairs, P.repo("config", "arm2_pairs.tsv"))
    P.log("arm 2: %d ScRAP strain names considered, %d pairs" % (len(cand), len(pairs)))

    # Extract the short-read assemblies of the pairs, in one pass over the archive.
    outdir = P.data_dir("assemblies", "peter2018")
    os.makedirs(outdir, exist_ok=True)
    need = [m for c, m in zip(pairs["code"], pairs["short_read_member"])
            if not os.path.exists(os.path.join(outdir, c + ".fna.gz"))]
    if need:
        work = P.data_dir("work", "peter_extract")
        os.makedirs(work, exist_ok=True)
        lst = os.path.join(work, "members.txt")
        P.atomic_write_text(lst, "\n".join(need) + "\n")
        subprocess.run(["tar", "-xzf", archive, "-C", work, "-T", lst], check=True)
        for c, m in zip(pairs["code"], pairs["short_read_member"]):
            src = os.path.join(work, m)
            if os.path.exists(src):
                subprocess.run("pigz -c '%s' > '%s.part' && mv '%s.part' '%s'" % (src, os.path.join(outdir, c + ".fna.gz"),
                               os.path.join(outdir, c + ".fna.gz"), os.path.join(outdir, c + ".fna.gz")),
                               shell=True, check=True)
                os.remove(src)
        subprocess.run(["rm", "-rf", work], check=True)
    P.log("short-read assemblies in %s" % P.rel(outdir))

    # Arm 3 runs: Peter et al. 2018 reads, linked through Table S17.
    runs = ena_runs()
    runs["stem"] = runs["submitted_ftp"].fillna("").str.split(";").map(
        lambda fs: sorted({submitted_stem(f) for f in fs if f}))
    s17i = s17.set_index("isolate")
    rows = []
    for pr in pairs.itertuples():
        code3 = pr.code.replace("SACE_", "")
        base = {"code": pr.code, "isolate_name": pr.isolate_name, "long_read_accession": pr.long_read_accession,
                "long_read_assembly_type": pr.long_read_assembly_type, "long_read_length": pr.long_read_length,
                "peter_zygosity": pr.peter_zygosity, "peter_aneuploidies": pr.peter_aneuploidies,
                "order_key": hashlib.sha256(("arm3|" + pr.code).encode()).hexdigest()}
        reads_file = s17i.loc[code3, "reads_file"] if code3 in s17i.index else None
        why = []
        if not str(pr.long_read_assembly_type).startswith("haploid/collapsed"):
            why.append("long-read assembly is not haploid or collapsed")
        if str(pr.peter_zygosity).strip().lower() != "homozygous":
            why.append("not homozygous in Peter et al. 2018")
        if str(pr.peter_aneuploidies).strip().lower() != "euploid":
            why.append("not euploid in Peter et al. 2018")
        linked = runs[runs["stem"].map(lambda s: isinstance(reads_file, str) and reads_file in s)] if reads_file else runs.iloc[0:0]
        linked = linked[(linked["library_layout"] == "PAIRED")]
        if not len(linked):
            why.append("no paired-end run in %s linked through Table S17" % ENA_STUDY)
            rows.append(dict(base, arm3_eligible="no", reason="; ".join(why)))
            continue
        r = linked.assign(bc=linked["base_count"].astype(float)).sort_values("bc", ascending=False).iloc[0]
        fq = (r["fastq_ftp"] or "").split(";")
        md5 = (r["fastq_md5"] or "").split(";")
        if len(fq) != 2:
            why.append("run does not have exactly two fastq files")
        rows.append(dict(base, run_accession=r["run_accession"], sample_accession=r["sample_accession"],
                         library_name=r["library_name"], instrument_model=r["instrument_model"],
                         read_count=r["read_count"], base_count=r["base_count"],
                         fastq_1=fq[0] if fq else "", fastq_2=fq[1] if len(fq) > 1 else "",
                         md5_1=md5[0] if md5 else "", md5_2=md5[1] if len(md5) > 1 else "", fastq_bytes=r["fastq_bytes"],
                         link_evidence="Table S17 reads file %s equals the submitted file names of %s"
                                       % (reads_file, r["run_accession"]),
                         depth_available=round(float(r["base_count"]) / float(pr.long_read_length), 1),
                         arm3_eligible="no" if why else "yes", reason="; ".join(why)))
    a3 = pd.DataFrame(rows).sort_values("order_key")
    P.atomic_write_tsv(a3, P.repo("config", "arm3_runs.tsv"))
    P.log("arm 3: %d strains with runs listed, %d eligible" % (len(a3), (a3["arm3_eligible"] == "yes").sum()))
    P.mark_done("build_pairs")


if __name__ == "__main__":
    main()
