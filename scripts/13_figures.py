"""Figures, one or more per question, from the tables in results/.

Colours are the first three slots of a palette checked for colour-blind
separation (blue, orange, aqua; worst all-pairs CVD distance 9.2); per-genome
or per-strain curves are thin grey lines behind the pooled one. Every panel
states its n and unit of replication. One y-axis per panel. Each figure is
written as PDF and PNG with no creation date in its metadata.

Usage:
    python scripts/13_figures.py
"""
import os
import sys
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import padlib as P  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GREY, LIGHT, INK, INK2, GRID = "#a3a29d", "#d9d8d3", "#0b0b0b", "#52514e", "#e6e5e1"
FIG = P.repo("figures")
SUPP = P.repo("figures", "supplementary")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y",
    "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-", "axes.axisbelow": True, "lines.linewidth": 1.5,
    "lines.solid_capstyle": "round", "pdf.fonttype": 42, "svg.hashsalt": "pcan", "figure.dpi": 100,
})


def read(path):
    p = P.repo("results", path)
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


def save(fig, name, supp=False):
    d = SUPP if supp else FIG
    os.makedirs(d, exist_ok=True)
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}),
                      ("png", {"Software": None})):
        fig.savefig(os.path.join(d, "%s.%s" % (name, ext)), dpi=300, bbox_inches="tight", metadata=meta)
    plt.close(fig)
    P.log("wrote %s" % P.rel(os.path.join(d, name + ".pdf")))


def end_label(ax, x, y, text, color=INK2, dy=0):
    ax.annotate(text, (x, y), xytext=(4, dy), textcoords="offset points", va="center", fontsize=7, color=color)


# ----------------------------------------------------------------------------
# Question 1, arm 0

def fig_arm0():
    causes = read("arm0/non_exact_cause_counts.tsv")
    expected = read("arm0/authors_expected_vs_observed.tsv")
    gate = read("arm0/gate.tsv")
    if causes is None or expected is None or gate is None:
        return
    fig = plt.figure(figsize=(7.2, 5.6), layout="constrained")
    top, bottom = fig.subfigures(2, 1, height_ratios=[1.1, 1], hspace=0.06)

    ax = top.add_subplot(1, 1, 1)
    c = causes.groupby("cause")["published_calls"].sum().sort_values()
    ax.barh(range(len(c)), c.values, height=0.6, color=BLUE)
    ax.set_yticks(range(len(c)))
    ax.set_yticklabels([s[0].upper() + s[1:].replace("_", " ") for s in c.index], fontsize=7)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    for i, v in enumerate(c.values):
        ax.text(v + 0.4, i, str(v), va="center", fontsize=7, color=INK2)
    ax.set_xlabel("published calls")
    n_sp = int((expected["published_calls"] > expected["exact"]).sum())
    ax.set_title("a  Why %d of %d published calls are not reproduced exactly (n = %d species with at least one; "
                 "unit: call)" % (c.sum(), int(gate.iloc[0]["n_calls"]), n_sp), loc="left")

    ax_b, ax_c = bottom.subplots(1, 2, width_ratios=[1.3, 1])
    ax = ax_b
    e = expected.dropna(subset=["authors_false_neg"]).copy()
    e["excess"] = (e["missing"] - e["authors_false_neg"]).astype(int)
    counts = e["excess"].value_counts().sort_index()
    ax.bar(counts.index, counts.values, width=0.6, color=BLUE)
    for x, v in counts.items():
        ax.text(x, v + 1.5, str(v), ha="center", fontsize=7, color=INK2)
    for _, r in e[e["excess"] >= 5].iterrows():
        ax.annotate(r["species"], (r["excess"], counts[r["excess"]]), xytext=(0, 14), textcoords="offset points",
                    ha="right", fontsize=6.5, color=INK2, style="italic")
    ticks = sorted(set([int(counts.index.min())] + list(range(0, int(counts.index.max()) + 1, 2))))
    ax.set_xticks(ticks)
    ax.set_ylim(0, counts.max() * 1.15)
    ax.set_xlabel("missing calls minus the false negatives\nthe authors' table predicts")
    ax.set_ylabel("species")
    ax.set_title("b  Misses against the authors' prediction\nn = %d species; unit: species" % len(e), loc="left")

    ax = ax_c
    g = gate.iloc[0]
    ax.errorbar([g["estimate"]], [0], xerr=[[g["estimate"] - g["ci_low"]], [g["ci_high"] - g["estimate"]]], fmt="o",
                color=BLUE, ecolor=BLUE, elinewidth=1.5, capsize=0, markersize=5)
    ax.axvline(0.95, color=ORANGE, linewidth=1)
    ax.annotate("pre-registered gate, 0.95", (0.95, 0.55), xytext=(4, 0), textcoords="offset points", fontsize=7,
                color=INK2)
    ax.annotate("%.3f\n(%.3f to %.3f)" % (g["estimate"], g["ci_low"], g["ci_high"]), (0.852, -0.45),
                ha="left", fontsize=7, color=INK2)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_xlim(0.85, 1.0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("fraction of published calls\nreproduced exactly")
    ax.set_title("c  The gate\nn = %d calls, %d species;\ninterval: bootstrap over species"
                 % (g["n_calls"], g["n_species"]), loc="left")
    save(fig, "fig1_reproduction")


# ----------------------------------------------------------------------------
# Question 2, arm 1a (simulation)

def fig_arm1a():
    per = read("arm1/fragmentation_replicates.tsv")
    stats0 = read("arm0/assembly_stats.tsv")
    stats2 = read("arm2/short_read_assembly_stats.tsv")
    stats3 = read("arm3/assemblies.tsv")
    if per is None:
        return
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), sharex=True)
    for col, model in enumerate(["uniform", "at_weighted"]):
        d = per[per["model"] == model]
        n_gen = d["accession"].nunique()
        n_rep = d.groupby(["accession", "target_n50"]).size().max()
        title = "uniform breakpoints" if model == "uniform" else "AT-weighted breakpoints"
        ax = axes[0, col]
        for acc, g in d.groupby("accession"):
            m = g.groupby("target_n50").agg(x=("realised_n50", "median"), y=("recall", "median")).sort_values("x")
            ax.plot(m["x"], m["y"], color=LIGHT, linewidth=0.8, zorder=1)
        lv = d.groupby("target_n50").agg(x=("realised_n50", "median"), med=("recall", "median"),
                                          q25=("recall", lambda s: s.quantile(0.25)),
                                          q75=("recall", lambda s: s.quantile(0.75)),
                                          geo=("geometric_expectation", "median"),
                                          syn=("synteny_checkable_fraction", "median")).sort_values("x")
        ax.fill_between(lv["x"], lv["q25"], lv["q75"], color=BLUE, alpha=0.12, linewidth=0, zorder=2)
        # The expectation is drawn wide and pale under the recall line: where
        # the two coincide, which is everywhere, only its edges show.
        ax.plot(lv["x"], lv["geo"], color=ORANGE, linewidth=5, alpha=0.45, solid_capstyle="round", zorder=3)
        ax.plot(lv["x"], lv["med"], color=BLUE, marker="o", markersize=3.5, zorder=4)
        ax.plot(lv["x"], lv["syn"], color=AQUA, marker="o", markersize=3.5, zorder=3)
        x_end = lv["x"].iloc[-1]
        end_label(ax, x_end, lv["med"].iloc[-1], "recall", dy=7)
        end_label(ax, x_end, lv["med"].iloc[-1], "geometric\nexpectation", dy=-12)
        ax.annotate("synteny-checkable", (lv["x"].iloc[3], lv["syn"].iloc[3]), xytext=(6, -2),
                    textcoords="offset points", fontsize=7, color=INK2)
        ax.annotate("%.2f" % lv["med"].iloc[0], (lv["x"].iloc[0], lv["med"].iloc[0]), xytext=(0, -12),
                    textcoords="offset points", ha="center", fontsize=7, color=INK2)
        if stats0 is not None:
            ax.plot(stats0["contig_n50"], np.full(len(stats0), -0.05), "|", color=GREY, markersize=6,
                    transform=ax.get_xaxis_transform(), clip_on=False)
        if stats2 is not None:
            ax.plot(stats2["contig_n50"], np.full(len(stats2), -0.10), "|", color=INK2, markersize=6,
                    transform=ax.get_xaxis_transform(), clip_on=False)
        if stats3 is not None and "contig_n50" in stats3:
            s3 = stats3.dropna(subset=["contig_n50"])
            ax.plot(s3["contig_n50"], np.full(len(s3), -0.15), "|", color=AQUA, markersize=6,
                    transform=ax.get_xaxis_transform(), clip_on=False)
        ax.set_xscale("log")
        ax.set_xlim(3000, 2.5e6)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("fraction of unfragmented calls")
        ax.set_title("%s  Simulation, %s\nn = %d genomes x %d replicates per level" % ("ab"[col], title, n_gen, n_rep),
                     loc="left")

        ax = axes[1, col]
        dec = d.groupby("target_n50").agg(x=("realised_n50", "median"), truth=("n_truth", "sum"),
                                           geo=("lost_geometric", "sum"), pipe=("lost_pipeline", "sum"),
                                           gained=("gained", "sum")).sort_values("x")
        # Labelled at the 5 kb end, where the three lines are apart; at the
        # other end all three meet at zero.
        for key, color, lab, dy in (("geo", ORANGE, "lost, window cut", 0), ("pipe", BLUE, "lost, window intact", 8),
                                    ("gained", AQUA, "gained", 8)):
            y = dec[key] / dec["truth"] * 100
            ax.plot(dec["x"], y, color=color, marker="o", markersize=3.5)
            ax.annotate(lab, (dec["x"].iloc[0], y.iloc[0]), xytext=(6, dy), textcoords="offset points",
                        fontsize=7, color=INK2, va="center")
        ax.set_xscale("log")
        ax.set_xlim(3000, 2.5e6)
        ax.set_xlabel("realised contig N50 (bp)")
        ax.set_ylabel("calls per 100 unfragmented calls")
        ax.set_title("%s  Simulation, losses and gains, %s\nn = %d genomes; calls pooled over replicates"
                     % ("cd"[col], title, n_gen), loc="left")
    ticks = ["the %d published arm 0 assemblies (light grey)" % (len(stats0) if stats0 is not None else 0)]
    if stats2 is not None:
        ticks.append("the %d published short-read assemblies of arm 2 (dark grey)" % len(stats2))
    if stats3 is not None and "contig_n50" in stats3:
        ticks.append("the %d arm 3 assemblies (aqua)" % len(stats3.dropna(subset=["contig_n50"])))
    note = textwrap.fill("Simulation. Recall is relative to PCAn's call on the unfragmented assembly. Lines: median "
                         "over all replicates; band: interquartile range; grey lines: each genome's median. Ticks "
                         "under the top axes, contig N50 of: %s." % "; ".join(ticks), 150)
    fig.text(0.01, -0.02, note, fontsize=6.8, color=INK2, ha="left", va="top")
    fig.tight_layout()
    save(fig, "fig2_fragmentation")

    # Supplementary: one panel per genome.
    gens = sorted(per["accession"].unique())
    fig, axes = plt.subplots(2, (len(gens) + 1) // 2, figsize=(7.2, 3.4), sharex=True, sharey=True)
    for ax, acc in zip(axes.flat, gens):
        for model, color in (("uniform", BLUE), ("at_weighted", ORANGE)):
            g = per[(per["accession"] == acc) & (per["model"] == model)]
            m = g.groupby("target_n50").agg(x=("realised_n50", "median"), y=("recall", "median")).sort_values("x")
            ax.plot(m["x"], m["y"], color=color, marker="o", markersize=2.5, linewidth=1)
        ax.set_xscale("log")
        ax.set_title(acc, fontsize=7)
        ax.set_ylim(0, 1.05)
    axes.flat[0].legend(["uniform", "AT-weighted"], frameon=False, fontsize=6.5)
    fig.supxlabel("realised contig N50 (bp)", fontsize=8)
    fig.supylabel("recall relative to the unfragmented call", fontsize=8)
    fig.suptitle("Simulation: recall per genome, median of 10 replicates per level; unit: genome", fontsize=8.5)
    fig.tight_layout()
    save(fig, "figS1_fragmentation_per_genome", supp=True)


# ----------------------------------------------------------------------------
# Question 3, arm 1b (simulation)

def fig_arm1b():
    env = read("arm1/detection_envelope.tsv")
    prog = read("arm1/planted_progressive.tsv")
    if env is None or prog is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    ax = axes[0]
    pooled = env[env["accession"] == "pooled"].sort_values("size")
    genomes = env[env["accession"] != "pooled"]
    for acc, g in genomes.groupby("accession"):
        g = g.sort_values("size")
        ax.plot(g["size"], g["correct_length_fraction"], color=LIGHT, linewidth=0.8)
    for key, lo, hi, color, lab in (("called_fraction", "called_ci_low", "called_ci_high", BLUE, "called"),
                                    ("correct_length_fraction", "correct_length_ci_low", "correct_length_ci_high",
                                     ORANGE, "called at the planted length")):
        ax.errorbar(pooled["size"], pooled[key], yerr=[pooled[key] - pooled[lo], pooled[hi] - pooled[key]],
                    color=color, marker="o", markersize=3.5, capsize=0, elinewidth=1, label=lab)
    ax.axvline(0, color=GRID, linewidth=0.8)
    ax.set_xticks(sorted(pooled["size"].unique()))
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("planted CDEII change (bp)")
    ax.set_ylabel("fraction of edited centromeres")
    ax.legend(frameon=False, loc="lower left")
    ax.set_title("a  Simulation: detection envelope\nn = %d edits in %d genomes;\ninterval: bootstrap over genomes; "
                 "grey: each genome" % (pooled["edits"].sum(), genomes["accession"].nunique()), loc="left")

    ax = axes[1]
    ok = prog[prog["status"] == "ok"]
    for key, color, lab, dx in (("total_calls", BLUE, "calls in the genome", -0.25),
                                ("edited_called", ORANGE, "edited centromeres called", 0),
                                ("edited_correct_length", AQUA, "called at +10 bp", 0.25)):
        ax.scatter(ok["k"] + dx, ok[key], s=10, color=color, edgecolor="white", linewidth=0.5, zorder=3, label=lab)
        med = ok.groupby("k")[key].median()
        ax.plot(med.index + dx, med.values, color=color, linewidth=1)
    ax.plot([0, 16], [0, 16], color=LIGHT, linewidth=0.8, zorder=1)
    ax.set_xticks(sorted(ok["k"].unique()))
    ax.set_ylim(0, 17.5)
    ax.set_xlabel("centromeres carrying a +10 bp insertion (S288C)")
    ax.set_ylabel("count")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("b  Simulation: a transition part of the way\nn = %d replicates per level; unit: replicate"
                 % ok.groupby("k").size().max(), loc="left")
    fig.tight_layout()
    save(fig, "fig3_planted_variants")


# ----------------------------------------------------------------------------
# Question 4, arm 2

STATUS_LABEL = {"intact_called": "intact, called", "intact_uncalled": "intact, not called", "split": "split",
                "n_run": "N run", "absent": "absent"}


def holm_p(source, claim_start):
    """The Holm-adjusted p-value of a registered test, from the reporting
    summary that scripts/12_analyse.py writes."""
    rep = read("reporting_summary.tsv")
    if rep is None or "p_holm" not in rep:
        return np.nan
    r = rep[(rep["source"] == source) & rep["claim"].astype(str).str.startswith(claim_start)]
    return float(r["p_holm"].iloc[0]) if len(r) else np.nan


def p_text(p):
    if p != p:
        return "see reporting summary"
    return "< 1e-16" if p < 1e-16 else "%.2g" % p


def fig_arm2():
    counts = read("arm2/call_counts.tsv")
    conf = read("arm2/confirmatory.tsv")
    cens = read("arm2/centromere_status.tsv")
    nulls = read("arm2/null_status.tsv")
    if counts is None or conf is None or cens is None or nulls is None:
        return
    n_str = counts["strain"].nunique()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))

    ax = axes[0, 0]
    diff = (counts["short_calls"] - counts["long_calls"]).value_counts().sort_index()
    ax.bar(diff.index, diff.values, width=0.6, color=BLUE)
    for x, v in diff.items():
        ax.text(x, v + 0.5, str(v), ha="center", fontsize=7, color=INK2)
    ax.set_xticks(range(int(diff.index.min()), int(diff.index.max()) + 1, 2))
    ax.set_xlabel("calls on the short-read assembly\nminus calls on the long-read assembly")
    ax.set_ylabel("strains")
    c1 = conf[conf["outcome"] == "C1"].iloc[0]
    ax.set_title("a  Call counts, %d of %d strains equal\n(%.2f, %.2f to %.2f); unit: strain"
                 % (c1["k"], c1["n_strains"], c1["estimate"], c1["ci_low"], c1["ci_high"]), loc="left")

    ax = axes[0, 1]
    c2 = conf[conf["outcome"] == "C2"].copy()
    c2["status"] = c2["measure"].str.replace("fraction ", "", regex=False)
    c2 = c2.set_index("status").reindex(list(STATUS_LABEL))
    y = np.arange(len(c2))[::-1]
    ax.errorbar(c2["estimate"], y, xerr=[c2["estimate"] - c2["ci_low"], c2["ci_high"] - c2["estimate"]], fmt="o",
                color=BLUE, ecolor=BLUE, capsize=0, markersize=4)
    for yy, (st, r) in zip(y, c2.iterrows()):
        ax.annotate("%d" % r["k"], (r["ci_high"], yy), xytext=(5, 0), textcoords="offset points", va="center",
                    fontsize=7, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([STATUS_LABEL[s] for s in c2.index])
    ax.set_xlim(0, 1)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("fraction of long-read centromeres")
    ax.set_title("b  Where each long-read centromere is\nn = %d centromeres in %d strains; interval over strains"
                 % (int(c2["n_centromeres"].iloc[0]), n_str), loc="left")

    ax = axes[1, 0]
    d = cens.loc[cens["status"] == "intact_called", "cdeii_difference"].dropna().astype(int)
    vc = d.value_counts().sort_index()
    ax.bar(vc.index, vc.values, width=0.7, color=BLUE)
    for x, v in vc.items():
        ax.text(x, v * 1.15, str(v), ha="center", fontsize=6.5, color=INK2)
    # Log scale: one bar (no difference) holds nearly all centromeres.
    ax.set_yscale("log")
    ax.set_ylim(0.7, vc.max() * 3)
    lo_x, hi_x = int(vc.index.min()), int(vc.index.max())
    ax.set_xticks(range(lo_x - lo_x % 2, hi_x + 1, 2))
    ax.set_xlabel("CDEII length, short-read call\nminus long-read call (bp)")
    ax.set_ylabel("centromeres (log scale)")
    c3 = conf[conf["outcome"] == "C3"].iloc[0]
    ax.set_title("c  CDEII length, %d of %d identical\nunit: centromere, %d strains" % (c3["k"], c3["n_centromeres"],
                                                                                    c3["n_strains"]), loc="left")

    ax = axes[1, 1]
    lift = cens[cens["status"] != "not_liftable"].copy()
    lift["broken"] = lift["status"].isin(["split", "n_run", "absent"])
    nulls = nulls.copy()
    nulls["broken"] = nulls["status"].isin(["split", "n_run", "absent"])
    per = pd.DataFrame({"cen": lift.groupby("strain")["broken"].mean(), "null": nulls.groupby("strain")["broken"].mean()})
    for _, r in per.iterrows():
        ax.plot([0, 1], [r["null"], r["cen"]], color=LIGHT, linewidth=0.7)
    ax.scatter(np.zeros(len(per)), per["null"], s=8, color=GREY, zorder=3)
    ax.scatter(np.ones(len(per)), per["cen"], s=8, color=GREY, zorder=3)
    c4 = conf[(conf["outcome"] == "C4") & conf["measure"].str.startswith("break-rate")].iloc[0]
    pooled = [c4["k_null"] / c4["n_null"], c4["k"] / c4["n_centromeres"]]
    ax.plot([0, 1], pooled, color=BLUE, marker="o", markersize=5, linewidth=1.5, zorder=4)
    for x, v in zip([0, 1], pooled):
        ax.annotate("%.3f" % v, (x, v), xytext=(-8 if x == 0 else 8, 0), textcoords="offset points",
                    ha="right" if x == 0 else "left", va="center", fontsize=7, color=INK2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["AT-matched\nnull windows", "centromeres"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("fraction broken in the short-read assembly")
    ax.set_title("d  Breaks: ratio %.2f (%.2f to %.2f), Holm p %s\nn = %d strains; grey: each strain; blue: pooled"
                 % (c4["estimate"], c4["ci_low"], c4["ci_high"],
                    p_text(holm_p("results/arm2/confirmatory.tsv", "C4: break-rate ratio")), c4["n_strains"]),
                 loc="left")
    fig.tight_layout()
    save(fig, "fig4_published_pairs")
    fig_arm2_zygosity(counts, lift, nulls)


def fig_arm2_zygosity(counts, lift, nulls):
    """Exploratory: the arm 2 measures split by Peter et al.'s zygosity call.
    Excluding heterozygous strains is the registered sensitivity analysis;
    this figure shows why it matters."""
    z = counts.set_index("strain")["zygosity"].str.lower()
    groups = [("homozygous", BLUE), ("heterozygous", ORANGE)]
    per = pd.DataFrame({"broken": lift.groupby("strain")["broken"].mean(),
                        "null": nulls.groupby("strain")["broken"].mean(),
                        "diff": (counts.set_index("strain")["short_calls"] - counts.set_index("strain")["long_calls"])})
    per["zygosity"] = z.reindex(per.index)
    rng = np.random.default_rng(20261009)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for ax, key, ylab, letter, title in (
            (axes[0], "broken", "fraction of centromeres broken", "a", "Centromeres broken per strain"),
            (axes[1], "diff", "short-read minus long-read calls", "b", "Call count difference per strain")):
        labels = []
        for i, (g, col) in enumerate(groups):
            v = per.loc[per["zygosity"] == g, key].dropna()
            x = i + rng.uniform(-0.12, 0.12, len(v))
            ax.scatter(x, v, s=10, color=col, edgecolor="white", linewidth=0.4, zorder=3)
            ax.plot([i - 0.25, i + 0.25], [v.median()] * 2, color=INK, linewidth=1.5, zorder=4)
            labels.append("%s\nn = %d strains" % (g, len(v)))
            if key == "broken":
                nv = per.loc[per["zygosity"] == g, "null"].median()
                ax.plot([i - 0.25, i + 0.25], [nv] * 2, color=GREY, linewidth=1.2, linestyle=(0, (2, 2)), zorder=2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(labels)
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylabel(ylab)
        ax.set_title("%s  %s\nunit: strain; black: median%s" % (
            letter, title, "\ngrey dashes: median for the null windows" if key == "broken" else ""), loc="left")
    axes[0].set_ylim(-0.08, 1)
    fig.text(0.01, -0.04, "Exploratory. Zygosity as Peter et al. 2018 (Table S1) report it. Arm 2, the published "
             "short-read assembly of each strain against its long-read assembly.", fontsize=6.8, color=INK2,
             ha="left", va="top")
    fig.tight_layout()
    save(fig, "figS2_zygosity", supp=True)


# ----------------------------------------------------------------------------
# Question 5, arm 3

DEPTH_ORDER = ["5", "10", "20", "40", "full"]
ASSEMBLERS = {"spades": (BLUE, "SPAdes"), "megahit": (ORANGE, "MEGAHIT")}


def depth_axis(ax):
    """Realised depth on a log axis, ticked at the target depths and the 80x cap."""
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([5, 10, 20, 40, 80]))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(3.5, 115)
    ax.set_xlabel("realised depth (x), median over strains")


def by_depth(conf, outcome, assembler, measure=None):
    d = conf[(conf["outcome"] == outcome) & (conf["assembler"] == assembler)
             & (conf["seed"].astype(str) == "11")].copy()
    if measure:
        d = d[d["measure"].str.startswith(measure)]
    d["target_depth"] = d["target_depth"].astype(str)
    return d.set_index("target_depth").reindex([o for o in DEPTH_ORDER if o in set(d["target_depth"])])


def fig_arm3():
    conf = read("arm3/confirmatory.tsv")
    asm = read("arm3/assemblies.tsv")
    if conf is None or asm is None or conf.empty:
        return
    asm["target_depth"] = asm["target_depth"].astype(str)
    ok = asm[asm["status"] == "ok"]
    depth = ok[ok["seed"].astype(str) == "11"].groupby(["assembler", "target_depth"])["realised_depth_raw"].median()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0))
    for ax, outcome, ylab, letter, title in (
            (axes[0, 0], "D1", "fraction of long-read calls intact and called", "a", "Recall relative to the long-read calls"),
            (axes[0, 1], "D2", "fraction with the long-read CDEII length", "b", "CDEII agreement among intact calls")):
        n = 0
        for i, (a, (col, lab)) in enumerate(ASSEMBLERS.items()):
            d = by_depth(conf, outcome, a)
            if d.empty:
                continue
            n = max(n, int(d["n_strains"].max()))
            x = [depth.get((a, t), np.nan) for t in d.index]
            ax.errorbar(x, d["estimate"], yerr=[d["estimate"] - d["ci_low"], d["ci_high"] - d["estimate"]],
                        color=col, marker="o", markersize=3.5, capsize=0, elinewidth=1)
            # the two assemblers can end at the same value; their labels are
            # set apart vertically so that both stay readable
            end_label(ax, x[-1], d["estimate"].iloc[-1], lab, dy=7 if i == 0 else -7)
            ten = conf[(conf["outcome"] == outcome) & (conf["assembler"] == a)
                       & (conf["target_depth"].astype(str) == "10")]
            if len(ten) > 1:
                ax.plot([depth.get((a, "10"), np.nan)] * len(ten), ten["estimate"], "_", color=col, markersize=9,
                        markeredgewidth=1.5)
        depth_axis(ax)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(ylab)
        ax.set_title("%s  %s\nn = %d strains; mean with bootstrap interval\nover strains; dashes at 10x: three seeds"
                     % (letter, title, n), loc="left")
    ax = axes[1, 0]
    for a, (col, lab) in ASSEMBLERS.items():
        d = ok[ok["assembler"] == a]
        ax.scatter(d["realised_depth_raw"], d["contig_n50"], s=9, color=col, edgecolor="white", linewidth=0.4,
                   label=lab, zorder=3)
    depth_axis(ax)
    ax.set_xlabel("realised depth (x)")
    ax.set_yscale("log")
    ax.set_ylabel("contig N50 (bp)")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("c  Contiguity, on the axis of arm 1\nn = %d assemblies; unit: assembly" % len(ok), loc="left")
    ax = axes[1, 1]
    for i, (a, (col, lab)) in enumerate(ASSEMBLERS.items()):
        d = by_depth(conf, "D4", a, "break-rate")
        if d.empty:
            continue
        x = [depth.get((a, t), np.nan) for t in d.index]
        ax.errorbar(x, d["estimate"], yerr=[d["estimate"] - d["ci_low"], d["ci_high"] - d["estimate"]], color=col,
                    marker="o", markersize=3.5, capsize=0, elinewidth=1)
        end_label(ax, x[-1], d["estimate"].iloc[-1], lab, dy=7 if i == 0 else -7)
    ax.axhline(1, color=GREY, linewidth=0.8)
    ax.annotate("no difference", (4, 1), xytext=(0, 3), textcoords="offset points", fontsize=7, color=INK2)
    depth_axis(ax)
    ax.set_ylabel("break-rate ratio, centromere over null")
    ax.set_title("d  Breaks at centromeres against AT-matched\nwindows, seed 11; bootstrap interval over strains",
                 loc="left")
    fig.tight_layout()
    save(fig, "fig5_read_depth")
    fig_arm3_per_strain(asm)


def fig_arm3_per_strain(asm):
    """Supplementary: recall per strain against realised depth, seed 11."""
    cen = read("arm3/centromere_status.tsv")
    if cen is None or cen.empty:
        return
    cen = cen[(cen["status"] != "not_liftable") & (cen["seed"].astype(str) == "11")].copy()
    cen["target_depth"] = cen["target_depth"].astype(str)
    rec = cen.groupby(["assembler", "strain", "target_depth"])["status"].apply(
        lambda s: (s == "intact_called").mean()).rename("recall").reset_index()
    rd = asm[asm["seed"].astype(str) == "11"][["assembler", "strain", "target_depth", "realised_depth_raw"]]
    rec = rec.merge(rd, on=["assembler", "strain", "target_depth"], how="left")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for ax, (a, (col, lab)) in zip(axes, ASSEMBLERS.items()):
        d = rec[rec["assembler"] == a]
        for _, g in d.groupby("strain"):
            g = g.sort_values("realised_depth_raw")
            ax.plot(g["realised_depth_raw"], g["recall"], color=LIGHT, linewidth=0.8, zorder=1)
        m = d.groupby("target_depth").agg(x=("realised_depth_raw", "median"), y=("recall", "mean")).reindex(
            [o for o in DEPTH_ORDER if o in set(d["target_depth"])])
        ax.plot(m["x"], m["y"], color=col, marker="o", markersize=3.5, zorder=3)
        depth_axis(ax)
        ax.set_xlabel("realised depth (x)")
        ax.set_ylim(0, 1.05)
        ax.set_title("%s  %s\nn = %d strains; grey: each strain; colour: mean"
                     % ("ab"[list(ASSEMBLERS).index(a)], lab, d["strain"].nunique()), loc="left")
    axes[0].set_ylabel("fraction of long-read calls\nintact and called")
    fig.tight_layout()
    save(fig, "figS3_read_depth_per_strain", supp=True)


# ----------------------------------------------------------------------------
# Across arms: where real assemblies break

def fig_breakpoints():
    c2 = read("arm2/confirmatory.tsv")
    c3 = read("arm3/confirmatory.tsv")
    bs = read("arm2/breakpoint_summary.tsv")
    if c2 is None or bs is None:
        return
    rows = []
    r = c2[c2["outcome"] == "C5"]
    if len(r):
        r = r.iloc[0]
        rows.append(("published short-read\nassemblies (arm 2)", r))
    if c3 is not None:
        for a in ("spades", "megahit"):
            r3 = c3[(c3["outcome"] == "D5") & (c3["assembler"] == a)]
            if len(r3):
                rows.append(("%s, full depth\n(arm 3)" % ("SPAdes" if a == "spades" else "MEGAHIT"), r3.iloc[0]))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1, 1.2]})
    ax = axes[0]
    d = (bs["mean_at_at_breakpoints"] - bs["genome_at"]).dropna() * 100
    bins = np.arange(np.floor(d.min() * 4) / 4, np.ceil(d.max() * 4) / 4 + 0.25, 0.25)
    hc, _, _ = ax.hist(d, bins=bins, color=BLUE, edgecolor="white", linewidth=0.6)
    ax.set_ylim(0, hc.max() * 1.2)
    ax.axvline(0, color=GREY, linewidth=0.8)
    ax.annotate("no AT preference", (0, ax.get_ylim()[1]), xytext=(3, -8), textcoords="offset points", ha="left",
                fontsize=7, color=INK2)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_xlabel("AT of windows holding a break minus genome AT\n(percentage points, 500 bp windows)")
    ax.set_ylabel("strains")
    ax.set_title("a  Published short-read assemblies (arm 2)\nn = %d strains, %d above zero; unit: strain"
                 % (len(d), int((d > 0).sum())), loc="left")
    ax = axes[1]
    y = np.arange(len(rows))[::-1]
    for yy, (lab, r) in zip(y, rows):
        ax.errorbar([r["estimate"]], [yy], xerr=[[r["estimate"] - r["ci_low"]], [r["ci_high"] - r["estimate"]]],
                    fmt="o", color=BLUE, capsize=0, markersize=4)
        ax.annotate("%.2f (%.2f to %.2f)\ndecision: %s" % (r["estimate"], r["ci_low"], r["ci_high"], r["decision"]),
                    (r["estimate"], yy), xytext=(0, -7), textcoords="offset points", ha="center", va="top",
                    fontsize=7, color=INK2)
    g1 = np.log(10) / 0.30
    ax.axvline(0, color=GREY, linewidth=0.8)
    ax.axvline(g1, color=ORANGE, linewidth=1)
    ax.annotate("uniform\nmodel", (0, len(rows) - 0.45), xytext=(4, 0), textcoords="offset points", fontsize=7,
                color=INK2, va="top")
    ax.annotate("arm 1 AT-\nweighted model", (g1, len(rows) - 0.45), xytext=(-4, 0), textcoords="offset points",
                fontsize=7, color=INK2, ha="right", va="top")
    ax.set_xlim(-0.8, max(g1, max(r["ci_high"] for _, r in rows)) + 2.5)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("maximum-likelihood gamma")
    ax.set_title("b  Which breakpoint model fits\ninterval: bootstrap over strains", loc="left")
    fig.tight_layout()
    save(fig, "fig6_breakpoint_model")


def main():
    fig_arm0()
    fig_arm1a()
    fig_arm1b()
    fig_arm2()
    fig_arm3()
    fig_breakpoints()


if __name__ == "__main__":
    main()
