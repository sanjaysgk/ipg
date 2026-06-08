#!/usr/bin/env python3
'''Self-contained cryptic-discovery HTML report.

Summarises the integrated peptide table for the cryptic subset: class
breakdown, cross-engine corroboration, PepQuery/spectral confidence, gffcompare
novelty class, RNA-seq origin, and before/after-rescore score shift split by
class. The modern concatenated-search analogue of the legacy db_compare Venn/
UpSet: one class column + cross-engine agreement replaces the two-DB Venn.
'''
from __future__ import annotations

import base64
import glob
import io
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SAMPLE = "${sample}"
PEPTIDES_TSV = "${peptides_tsv}"
CANONICAL_PREFIXES = [p for p in "${canonical_prefixes}".split(",") if p]
PROCESS_NAME = "${task.process}"

CLASS_COLORS = {"cryptic": "firebrick", "canonical": "steelblue"}


def fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def embed(fig, alt="plot") -> str:
    return f'<img alt="{alt}" src="data:image/png;base64,{fig_to_b64(fig)}" />'


def section(title, body) -> str:
    return f"<section><h2>{title}</h2>{body}</section>"


def classify(protein_list, prefixes) -> str:
    s = str(protein_list).strip("[]").replace("'", "").replace(", ", ";")
    prots = [p[4:] if p.startswith("rev_") else p for p in s.split(";") if p]
    if any(p.startswith(pre) for p in prots for pre in prefixes):
        return "canonical"
    return "cryptic"


df = pd.read_csv(PEPTIDES_TSV, sep="\\t")
if "class" not in df.columns:
    df["class"] = "canonical"
cryptic = df[df["class"] == "cryptic"]
n_cryptic, n_canon = len(cryptic), int((df["class"] == "canonical").sum())

parts = [f"<p><b>{n_canon}</b> canonical &middot; <b>{n_cryptic}</b> cryptic "
         f"peptides (total {len(df)}).</p>"]

# 1. Class breakdown.
fig, ax = plt.subplots(figsize=(4, 3.5))
ax.bar(["canonical", "cryptic"], [n_canon, n_cryptic],
       color=[CLASS_COLORS["canonical"], CLASS_COLORS["cryptic"]], edgecolor="black")
ax.set_ylabel("Peptides")
ax.set_title("Class breakdown")
parts.append(section("Class breakdown", embed(fig, "class")))

# 2. Cross-engine corroboration by class.
if "#engine" in df.columns:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    pivot = (df.groupby(["#engine", "class"]).size().unstack(fill_value=0))
    pivot.plot(kind="bar", ax=ax,
               color=[CLASS_COLORS.get(c, "gray") for c in pivot.columns],
               edgecolor="black")
    ax.set_xlabel("# engines supporting the peptide")
    ax.set_ylabel("Peptides")
    ax.set_title("Cross-engine corroboration")
    parts.append(section("Cross-engine corroboration (UpSet analogue)", embed(fig, "engine")))

# 3. Cryptic confidence: PepQuery + spectral.
conf_body = ""
for col, label in [("pepquery_status", "PepQuery2"), ("spectral_status", "MS2PIP/DeepLC")]:
    if col in cryptic.columns and not cryptic.empty:
        vc = cryptic[col].replace("", "n/a").value_counts()
        fig, ax = plt.subplots(figsize=(4.5, 3))
        ax.bar(vc.index.astype(str), vc.values, color="firebrick", edgecolor="black")
        ax.set_ylabel("Cryptic peptides")
        ax.set_title(f"{label} status")
        conf_body += embed(fig, col)
if conf_body:
    parts.append(section("Cryptic confidence gates", conf_body))

# 4. Novelty class (gffcompare) for cryptic.
if "novelty_class" in cryptic.columns and not cryptic.empty:
    vc = cryptic["novelty_class"].replace("", "n/a").value_counts()
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.bar(vc.index.astype(str), vc.values, color="darkorange", edgecolor="black")
    ax.set_xlabel("gffcompare class code")
    ax.set_ylabel("Cryptic peptides")
    ax.set_title("RNA-seq novelty class")
    parts.append(section("Cryptic origin — novelty class", embed(fig, "novelty")))

# 5. Before/after-rescore score shift, split by class (per-engine TSVs).
frames = []
for tsv in sorted(glob.glob("rescored/*/*.tsv") + glob.glob("rescored/*.tsv")):
    try:
        r = pd.read_csv(tsv, sep="\\t")
    except Exception:
        continue
    if not {"score", "provenance:before_rescoring_score", "is_decoy",
            "protein_list"}.issubset(r.columns):
        continue
    r = r[~r["is_decoy"].astype(bool)].copy()
    r["class"] = r["protein_list"].apply(lambda x: classify(x, CANONICAL_PREFIXES))
    frames.append(r[["class", "score", "provenance:before_rescoring_score"]])
if frames:
    allr = pd.concat(frames, ignore_index=True)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=False)
    for ax, cls in zip(axes, ["canonical", "cryptic"]):
        sub = allr[allr["class"] == cls]
        if not sub.empty:
            ax.hist(sub["provenance:before_rescoring_score"], bins=30, alpha=0.5,
                    label="before", color="gray")
            ax.hist(sub["score"], bins=30, alpha=0.5, label="after",
                    color=CLASS_COLORS.get(cls, "black"))
        ax.set_title(f"{cls} (n={len(sub)})")
        ax.set_xlabel("PSM score")
        ax.legend()
    fig.suptitle("Rescoring shift, before vs after, by class")
    parts.append(section("Rescoring gain by class", embed(fig, "rescore")))

# Cryptic discovery table (provenance + confidence).
table_cols = [c for c in ["peptide", "#engine", "peptide_qval_lowest",
              "pepquery_status", "spectral_status", "gene_name", "genomic_locus",
              "novelty_class", "FPKM", "TPM"] if c in cryptic.columns]
if not cryptic.empty and table_cols:
    parts.append(section("Cryptic discovery table",
                         cryptic[table_cols].to_html(index=False, na_rep="")))

# Summary TSV.
summary_rows = [("class:canonical", n_canon), ("class:cryptic", n_cryptic)]
for col in ["pepquery_status", "spectral_status", "novelty_class"]:
    if col in cryptic.columns:
        for k, v in cryptic[col].replace("", "n/a").value_counts().items():
            summary_rows.append((f"{col}:{k}", int(v)))
pd.DataFrame(summary_rows, columns=["metric", "count"]).to_csv(
    f"{SAMPLE}_cryptic_summary.tsv", sep="\\t", index=False)

html = (f"<html><head><meta charset='utf-8'><title>Cryptic report — {SAMPLE}</title>"
        "<style>body{font-family:sans-serif;margin:2em;max-width:1000px}"
        "section{margin:1.5em 0}img{max-width:100%}"
        "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:3px 6px}"
        "</style></head><body>"
        f"<h1>Cryptic peptide discovery report</h1><p>Sample: <b>{SAMPLE}</b></p>"
        + "".join(parts) + "</body></html>")
with open(f"{SAMPLE}_cryptic_report.html", "w") as fh:
    fh.write(html)
print(f"[cryptic_report] {SAMPLE}: {n_cryptic} cryptic / {len(df)} peptides",
      file=sys.stderr)

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    f.write(f"    pandas: {pd.__version__}\\n")
    f.write(f"    matplotlib: {matplotlib.__version__}\\n")
