"""Convert the downloaded supplementary tables to TSV, compare the publisher and
Figshare copies, and build config/species_arm0.tsv.

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
EFP_FIELDS = ["CDEIIImotif", "CDEImotif", "CDEIIIthreshold", "CDEIthreshold", "Works?", "FalseNeg", "FalsePos"]


def read_authors_table(path, known_names):
    """PCAn's ExpectedFalsePosFalseNegs.txt was saved with its tabs expanded to
    spaces, so a single space separates some columns as well as the words of a
    species name. Each line is index, genus, species name, then the seven
    fields above, the last two blank for one species. The species name is taken
    as the longest run of words, after the genus, that is a name in
    Supplementary Data 5; a line that matches none stops the run."""
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            name = None
            for k in range(min(5, len(tok) - 2), 1, -1):
                cand = " ".join(tok[2:2 + k])
                if cand in known_names:
                    name, rest = cand, tok[2 + k:]
                    break
            if name is None:
                raise SystemExit("cannot parse a species name from: %s" % line.strip())
            rest = rest + [""] * (len(EFP_FIELDS) - len(rest))
            row = {"SpeciesIndex": tok[0], "Genus": tok[1], "NewSpeciesName": name}
            row.update(dict(zip(EFP_FIELDS, rest)))
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    sd2 = read_sheet("supp_data2_publisher.xlsx")
    sd4 = read_sheet("supp_data4_publisher.xlsx")
    sd5 = read_sheet("supp_data5_publisher.xlsx")
    sd6_summary, sd6 = read_sd6("supp_data6_publisher.xlsx")
    for df, name in [(sd2, "supp_data2.tsv"), (sd4, "supp_data4.tsv"), (sd5, "supp_data5.tsv"),
                     (sd6, "supp_data6_assemblies.tsv"), (sd6_summary, "supp_data6_summary.tsv")]:
        P.atomic_write_tsv(df, os.path.join(T, name))

    # Publisher and Figshare copies. Supplementary Data 2 on Figshare is at
    # version 2 and differs in size from the publisher's file, so the two are
    # compared cell by cell and the publisher's copy, which is the version of
    # record, is the one used throughout.
    rows = []
    pairs = [("Supplementary Data 2", sd2, read_sheet("supp_data2_figshare_28225061v2.xlsx")),
             ("Supplementary Data 5", sd5, read_sheet("supp_data5_figshare_28225067v3.xlsx")),
             ("Supplementary Data 6", sd6, read_sd6("supp_data6_figshare_28225076v1.xlsx")[1])]
    for label, a, b in pairs:
        res = compare(a, b)
        if isinstance(res, tuple):
            summary, diff = res
            if summary["differing_cells"]:
                cols = [c for c in diff.columns if diff[c].any()]
                summary["differing_columns"] = ",".join(cols)
                if label == "Supplementary Data 2":
                    idx = diff.any(axis=1)
                    detail = pd.concat([normalise(a)[idx].add_suffix("_publisher"),
                                        normalise(b)[idx].add_suffix("_figshare")], axis=1)
                    P.atomic_write_tsv(detail, P.repo("results", "arm0", "source_differences_supp_data2.tsv"))
            else:
                summary["differing_columns"] = ""
        else:
            summary = res
            summary["differing_columns"] = "shapes differ"
        summary["table"] = label
        rows.append(summary)
    conc = pd.DataFrame(rows)[["table", "rows_a", "rows_b", "same_shape", "cells", "differing_cells", "differing_columns"]]
    conc = conc.rename(columns={"rows_a": "publisher_rows", "rows_b": "figshare_rows"})
    P.atomic_write_tsv(conc, P.repo("results", "arm0", "source_concordance.tsv"))

    # PCAn's own settings and the authors' per-species table.
    st = P.pcan_settings()
    genera = st["available_genus_list"]
    kaz = st["available_kazachstania_species"]
    mt = st["MotifsAndThresholds"]
    efp = read_authors_table(os.path.join(P.pcan_dir(), "ExpectedFalsePosFalseNegs.txt"),
                             set(sd5["Species_Name"].astype(str).str.strip()))
    if len(efp) != 138 or efp["NewSpeciesName"].duplicated().any():
        raise SystemExit("expected 138 distinct species in the authors' table, found %d" % len(efp))
    efp = efp.set_index("NewSpeciesName")

    def short(path, prefix):
        b = os.path.basename(path)
        return b[len(prefix):-len("_MEME.txt")] if b.startswith(prefix) and b.endswith("_MEME.txt") else b

    sd2_by_acc = sd2.groupby("Genbank_Assembly")
    sd2_acc_by_name = sd2.groupby("Species")["Genbank_Assembly"].first().to_dict()
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
            row["in_authors_table"] = "yes" if sp in efp.index else "no"
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
                           pcan_cdeiii_thresh="1e-%d" % t3, pcan_cdei_thresh="1e-%d" % t1)
            if sp in efp.index:
                e = efp.loc[sp]
                row.update(authors_cdeiii_motif=e["CDEIIImotif"], authors_cdei_motif=e["CDEImotif"],
                           authors_cdeiii_thresh=("1e-" + e["CDEIIIthreshold"]) if str(e["CDEIIIthreshold"]).isdigit() else e["CDEIIIthreshold"],
                           authors_cdei_thresh=("1e-" + e["CDEIthreshold"]) if str(e["CDEIthreshold"]).isdigit() else e["CDEIthreshold"],
                           authors_works=e["Works?"], authors_false_neg=e["FalseNeg"], authors_false_pos=e["FalsePos"])
                same = all(str(row.get("pcan_" + k, "")) == str(row.get("authors_" + k, ""))
                           for k in ("cdeiii_motif", "cdei_motif", "cdeiii_thresh", "cdei_thresh"))
                row["settings_match_authors_table"] = "yes" if same else "no"
            # Supplementary Data 2 is matched by assembly accession first,
            # because one species carries a different name there
            # (Vanderwaltozyma yanomamii appears as Vanderwaltozyma sp.
            # UFMG-CM-Y7006 under the same accession), and by name second,
            # because one assembly carries a different version there
            # (Jamesozyma spencerorum, GCA_003708825.3 in Data 5 and .2 in Data 2).
            if acc in sd2_by_acc.groups:
                calls = sd2_by_acc.get_group(acc)
                row["sd2_accession"] = acc
            elif sp in sd2_acc_by_name:
                calls = sd2[sd2["Species"] == sp]
                row["sd2_accession"] = sd2_acc_by_name[sp]
            else:
                calls = sd2.iloc[0:0]
                row["sd2_accession"] = ""
            row["sd2_name"] = calls["Species"].iloc[0] if len(calls) else ""
            row["sd2_calls"] = len(calls)
            row["sd2_calls_with_cdei"] = int((calls["CDEI"].astype(str).str.len() == 8).sum()) if len(calls) else 0
            if genus == "Naumovozyma":
                row.update(arm0_role="excluded",
                           reason="published calls came from an adapted single-motif version of PCAn "
                                  "(Extended Data Fig. 2) that the repository does not provide; PCAn v1.0 "
                                  "has no Naumovozyma entry in its motif table")
            elif key not in mt:
                row.update(arm0_role="excluded", reason="no entry in PCAn v1.0's motif table")
            elif len(calls) == 0:
                row.update(arm0_role="no_published_calls",
                           reason="no rows in Supplementary Data 2 and no motif in the authors' table; "
                                  "run and reported, outside the reproduction fraction")
            elif row["sd2_calls_with_cdei"] == 0:
                row.update(arm0_role="cdeiii_only",
                           reason="published calls have no CDEI (written ?) and no CDEII length; the released "
                                  "two-motif pipeline always reports a CDEI, so only CDEIII loci are compared, "
                                  "outside the reproduction fraction")
            else:
                row.update(arm0_role="primary", reason="")
        out.append(row)

    cols = ["sd5_row", "species", "accession", "group", "arm0_role", "reason", "pcan_genus",
            "pcan_kazachstania_species", "pcan_cdeiii_motif", "pcan_cdeiii_thresh", "pcan_cdei_motif",
            "pcan_cdei_thresh", "in_authors_table", "authors_cdeiii_motif", "authors_cdeiii_thresh",
            "authors_cdei_motif", "authors_cdei_thresh", "settings_match_authors_table", "authors_works",
            "authors_false_neg", "authors_false_pos", "sd2_name", "sd2_accession", "sd2_calls",
            "sd2_calls_with_cdei"]
    species = pd.DataFrame(out).reindex(columns=cols)
    P.atomic_write_tsv(species, P.repo("config", "species_arm0.tsv"))
    n = species["arm0_role"].value_counts().to_dict()
    P.log("config/species_arm0.tsv: %d rows, roles %s" % (len(species), n))
    P.log("Saccharomycetaceae: %d" % (species["group"] == "Saccharomycetaceae").sum())


if __name__ == "__main__":
    main()
