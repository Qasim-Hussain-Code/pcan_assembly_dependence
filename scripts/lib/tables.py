"""Convert the downloaded tables to TSV, compare the publisher and Figshare
copies, and build config/species_arm0.tsv.

Called by scripts/02_fetch_tables.sh. Every rule that decides a species' role in
arm 0 is written out below with its reason, and the reason is copied into the
species table so that no exclusion is silent.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import padlib as P  # noqa: E402

T = P.data_dir("tables")


def read_sheet(name, header_row=1):
    # Every supplementary sheet starts with an empty first row and an empty
    # first column.
    df = pd.read_excel(os.path.join(T, name), header=header_row)
    df = df.drop(columns=[c for c in df.columns if str(c).startswith("Unnamed: 0")])
    return df


def read_sd6(name):
    raw = pd.read_excel(os.path.join(T, name), header=None)
    raw = raw.drop(columns=[0])
    # Two blocks: a summary of the five sources, then one row per assembly
    # below a second header row reading Dataset, Assembly.
    second = raw.index[(raw[1] == "Dataset") & (raw[2] == "Assembly")][0]
    summary = raw.iloc[1:second].dropna(how="all")
    summary.columns = list(raw.iloc[0])
    summary = summary.dropna(axis=1, how="all")
    assemblies = raw.iloc[second + 1:, :2].dropna(how="all")
    assemblies.columns = ["Dataset", "Assembly"]
    return summary.reset_index(drop=True), assemblies.reset_index(drop=True)


def normalise(df):
    out = df.copy()
    for c in out.columns:
        out[c] = out[c].map(lambda v: "" if pd.isna(v) else str(v).strip())
    return out.reset_index(drop=True)


def compare(a, b):
    """Cell-by-cell agreement of two copies of one table."""
    a, b = normalise(a), normalise(b)
    same_shape = a.shape == b.shape and list(a.columns) == list(b.columns)
    if not same_shape:
        return {"same_shape": "no", "cells": "", "differing_cells": "", "rows_a": len(a), "rows_b": len(b)}
    diff = (a != b)
    return {"same_shape": "yes", "cells": int(diff.size), "differing_cells": int(diff.values.sum()),
            "rows_a": len(a), "rows_b": len(b)}, diff


OUTGROUPS = {"Wickerhamomyces anomalus", "Candida albicans", "Pichia kudriavzevii", "Yarrowia lipolytica"}
MUCOROMYCOTA_GENERA = {"Mucor", "Parasitella", "Ellisomyces", "Chaetocladium"}


def main():
    sd4 = read_sheet("supp_data4_publisher.xlsx")
    sd5 = read_sheet("supp_data5_publisher.xlsx")
    sd6_summary, sd6 = read_sd6("supp_data6_publisher.xlsx")
    for df, name in [(sd4, "supp_data4.tsv"), (sd5, "supp_data5.tsv"),
                     (sd6, "supp_data6_assemblies.tsv"), (sd6_summary, "supp_data6_summary.tsv")]:
        P.atomic_write_tsv(df, os.path.join(T, name))

    # Publisher and Figshare copies, compared cell by cell. The publisher's
    # copy is the one used throughout.
    rows = []
    pairs = [("Supplementary Data 5", sd5, read_sheet("supp_data5_figshare_28225067v3.xlsx")),
             ("Supplementary Data 6", sd6, read_sd6("supp_data6_figshare_28225076v1.xlsx")[1])]
    for label, a, b in pairs:
        res = compare(a, b)
        if isinstance(res, tuple):
            summary, diff = res
            summary["differing_columns"] = ",".join(c for c in diff.columns if diff[c].any())
        else:
            summary = res
            summary["differing_columns"] = "shapes differ"
        summary["table"] = label
        rows.append(summary)
    conc = pd.DataFrame(rows)[["table", "rows_a", "rows_b", "same_shape", "cells", "differing_cells", "differing_columns"]]
    conc = conc.rename(columns={"rows_a": "publisher_rows", "rows_b": "figshare_rows"})
    P.atomic_write_tsv(conc, P.repo("results", "arm0", "source_concordance.tsv"))

    # PCAn's own genus list and motif table decide which species run.
    st = P.pcan_settings()
    genera = st["available_genus_list"]
    mt = st["MotifsAndThresholds"]

    def short(path, prefix):
        b = os.path.basename(path)
        return b[len(prefix):-len("_MEME.txt")] if b.startswith(prefix) and b.endswith("_MEME.txt") else b

    out = []
    for i, r in sd5.reset_index(drop=True).iterrows():
        sp = str(r["Species_Name"]).strip()
        acc = str(r["Genbank_Assembly"]).strip()
        genus = sp.split()[0]
        row = {"sd5_row": i + 1, "species": sp, "accession": acc}
        if sp in OUTGROUPS:
            row.update(group="outgroup", arm0_role="excluded",
                       reason="outgroup, not Saccharomycetaceae; out of scope for arm 0")
        elif genus == "Hanseniaspora":
            row.update(group="Saccharomycodaceae", arm0_role="excluded",
                       reason="Saccharomycodaceae, not Saccharomycetaceae; out of scope for arm 0")
        elif genus in MUCOROMYCOTA_GENERA:
            row.update(group="Mucoromycota", arm0_role="excluded",
                       reason="Mucoromycota; predicted with PCAn_mucor, out of scope for arm 0")
        else:
            row["group"] = "Saccharomycetaceae"
            if genus not in genera:
                row.update(arm0_role="excluded", reason="genus not in PCAn v1.0's genus list")
                out.append(row)
                continue
            key = genus
            if genus == "Kazachstania":
                epithet = sp.split()[1]
                row["pcan_kazachstania_species"] = epithet
                key = "Kazachstania " + epithet
            row["pcan_genus"] = genus
            if key in mt:
                c3, c1, t3, t1 = mt[key]
                row.update(pcan_cdeiii_motif=short(c3, "CDEIII_"), pcan_cdei_motif=short(c1, "CDEI_"),
                           pcan_cdeiii_thresh="1e-%d" % t3, pcan_cdei_thresh="1e-%d" % t1,
                           arm0_role="run", reason="")
            else:
                row.update(arm0_role="excluded", reason="no entry in PCAn v1.0's motif table")
        out.append(row)

    cols = ["sd5_row", "species", "accession", "group", "arm0_role", "reason", "pcan_genus",
            "pcan_kazachstania_species", "pcan_cdeiii_motif", "pcan_cdeiii_thresh", "pcan_cdei_motif",
            "pcan_cdei_thresh"]
    species = pd.DataFrame(out).reindex(columns=cols)
    P.atomic_write_tsv(species, P.repo("config", "species_arm0.tsv"))
    n = species["arm0_role"].value_counts().to_dict()
    P.log("config/species_arm0.tsv: %d rows, roles %s" % (len(species), n))
    P.log("Saccharomycetaceae: %d" % (species["group"] == "Saccharomycetaceae").sum())


if __name__ == "__main__":
    main()
