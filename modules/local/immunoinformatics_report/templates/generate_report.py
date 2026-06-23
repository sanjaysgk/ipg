#!/usr/bin/env python3
'''Generate a self-contained immunoinformatics HTML report.

Consumes per-sample outputs from the IMMUNOINFORMATICS subworkflow and
produces a single HTML page with embedded PNG plots. Inputs are all
optional — missing tables degrade gracefully to "section skipped".
Sequence logos use logomaker when available.
'''
from __future__ import annotations

import argparse
import base64
import io
import sys
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

try:
    import logomaker as lm
    HAS_LOGOMAKER = True
except ImportError:
    HAS_LOGOMAKER = False


BIND_COLORS = {"SB": "firebrick", "WB": "dodgerblue",
               "NB": "gray", "NA": "gainsboro"}


def fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def embed(png_b64: str, alt: str = "plot") -> str:
    return f'<img alt="{alt}" src="data:image/png;base64,{png_b64}" />'


def section(title: str, body: str) -> str:
    return f'<section><h2>{title}</h2>{body}</section>'


def plot_length_histogram(peptides: pd.DataFrame) -> str:
    if "length" not in peptides.columns or peptides.empty:
        return "<p><em>No peptides to plot.</em></p>"
    lens = range(int(peptides["length"].min()), int(peptides["length"].max()) + 1)
    counts = peptides["length"].value_counts().reindex(lens, fill_value=0)
    immuno_mask = peptides.get("immuno", pd.Series([False] * len(peptides))).astype(bool)
    colors = ["firebrick" if immuno_mask[peptides["length"] == L].any() else "gray"
              for L in counts.index]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(counts.index.astype(str), counts.values,
           color=colors, edgecolor="black")
    ax.set_xlabel("Peptide length")
    ax.set_ylabel("Peptides")
    total = len(peptides)
    immuno_total = int(immuno_mask.sum())
    pct = (immuno_total / total * 100) if total else 0
    ax.set_title(f"Length histogram — {immuno_total}/{total} "
                 f"({pct:.1f}%) immunopeptides")
    return embed(fig_to_b64(fig), alt="length histogram")


def plot_psms_per_run(peptides: pd.DataFrame) -> str:
    run_cols = [c for c in peptides.columns if c.startswith("PSMs_run_")]
    if not run_cols:
        return "<p><em>No per-run PSM columns found.</em></p>"

    # Pull bind_level from peptides if present — else label everything NA.
    peptides = peptides.copy()
    peptides["bind_level"] = peptides.get("bind_level", "NA").fillna("NA")
    levels = ["SB", "WB", "NB", "NA"]
    psm_data = pd.DataFrame(0, index=run_cols, columns=levels)
    pep_data = pd.DataFrame(0, index=run_cols, columns=levels)
    for _, row in peptides.iterrows():
        lvl = row["bind_level"] if row["bind_level"] in levels else "NA"
        for col in run_cols:
            n = row.get(col, 0) or 0
            if n > 0:
                psm_data.loc[col, lvl] += n
                pep_data.loc[col, lvl] += 1

    psm_data.index = psm_data.index.str.replace("PSMs_run_", "")
    pep_data.index = pep_data.index.str.replace("PSMs_run_", "")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, data, ylabel in zip(axes,
                                (psm_data, pep_data),
                                ("PSMs", "Peptides")):
        data.plot(kind="bar", stacked=True, ax=ax,
                  color=[BIND_COLORS[c] for c in data.columns],
                  edgecolor="black")
        ax.set_xlabel("Run")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} per run")
        ax.legend(title="Binding level", fontsize="small")
    fig.tight_layout()
    return embed(fig_to_b64(fig), alt="psms per run")


def plot_netmhcpan_overview(best: pd.DataFrame, title: str) -> str:
    if best.empty or "bind_level" not in best.columns:
        return ""
    lengths = best["peptide"].str.len()
    best = best.assign(length=lengths)
    lens = range(int(best["length"].min()), int(best["length"].max()) + 1)
    counts = (best.groupby(["length", "bind_level"]).size()
              .unstack(fill_value=0).reindex(lens, fill_value=0))
    existing = [lvl for lvl in ["SB", "WB", "NB", "NA"] if lvl in counts.columns]
    counts = counts[existing]

    fig, ax = plt.subplots(figsize=(9, 5))
    counts.plot(kind="bar", stacked=True, ax=ax,
                color=[BIND_COLORS[c] for c in counts.columns],
                edgecolor="black")
    ax.set_xlabel("Peptide length")
    ax.set_ylabel("Peptides")
    ax.set_title(title)
    fig.tight_layout()
    return embed(fig_to_b64(fig), alt=title)


def plot_gibbs_logos(clusters: pd.DataFrame) -> str:
    if not HAS_LOGOMAKER or clusters.empty or "gibbs_core" not in clusters.columns:
        return "<p><em>logomaker unavailable or no cluster data.</em></p>"
    uniq = sorted(clusters["cluster"].dropna().astype(int).unique())
    if not uniq:
        return "<p><em>No Gibbs clusters found.</em></p>"

    fig, axes = plt.subplots(len(uniq), 1, figsize=(10, 3 * len(uniq)))
    if len(uniq) == 1:
        axes = [axes]
    for ax, cid in zip(axes, uniq):
        cores = clusters[clusters["cluster"] == cid]["gibbs_core"].astype(str).tolist()
        if len(cores) < 5:
            ax.set_title(f"Cluster {cid} (n={len(cores)}, too few to plot)")
            ax.axis("off")
            continue
        counts = lm.alignment_to_matrix(cores)
        info = lm.transform_matrix(counts, from_type="counts",
                                   to_type="information")
        lm.Logo(info, ax=ax, color_scheme="weblogo_protein")
        ax.set_title(f"Cluster {cid} (n={len(cores)})")
    fig.tight_layout()
    return embed(fig_to_b64(fig), alt="gibbs logos")


def quant_summary(quant: pd.DataFrame) -> str:
    if quant.empty:
        return "<p><em>No quantification data.</em></p>"
    n_pep = len(quant)
    run_cols = [c for c in quant.columns if "Intensity_" in c]
    n_runs = len(run_cols)
    return (f"<p>Quantified <strong>{n_pep}</strong> peptides across "
            f"<strong>{n_runs}</strong> run(s).</p>"
            + quant.head(20).to_html(index=False, classes="quant-preview"))


def read_tsv(path: Optional[Path]) -> pd.DataFrame:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\\t")
    except Exception as e:
        print(f"[report] failed to parse {path}: {e}", file=sys.stderr)
        return pd.DataFrame()


HTML_TEMPLATE = '''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>IPG immunoinformatics report — {sample}</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, sans-serif; max-width: 960px;
          margin: 2em auto; padding: 0 1em; color: #222; }}
  h1 {{ border-bottom: 2px solid #444; padding-bottom: .3em; }}
  h2 {{ margin-top: 2em; border-bottom: 1px solid #ccc; padding-bottom: .2em; }}
  img {{ max-width: 100%; height: auto; display: block; margin: 1em 0; }}
  table.quant-preview {{ border-collapse: collapse; font-size: .85em; }}
  table.quant-preview th, table.quant-preview td {{
    border: 1px solid #ddd; padding: 3px 8px; }}
  .meta {{ color: #666; font-size: .9em; }}
</style>
</head>
<body>
<h1>IPG immunoinformatics report — {sample}</h1>
<p class="meta">Sample: <code>{sample}</code> &middot; generated by sanjaysgk/ipg</p>
{summary}
{sections}
</body>
</html>
'''


# --- Nextflow template entry point ---------------------------------
# Optional inputs come in as 'NO_FILE' sentinel paths when the upstream
# module didn't run. _opt() converts those into None.
SAMPLE = "${meta.id}"
PEPTIDES = "${peptides}"
NETMHCPAN_RAW = "${netmhcpan}"
NETMHCIIPAN_RAW = "${netmhciipan}"
GIBBS_RAW = "${gibbs}"
FLASHLFQ_RAW = "${flashlfq}"
BLASTP_RAW = "${blastp}"
OUT = "${meta.id}_immunoinformatics_report.html"
PEPTIDE_LENGTH = "${pep_len}"
PROCESS_NAME = "${task.process}"


def _opt(p: str) -> Optional[Path]:
    return None if Path(p).name == "NO_FILE" else Path(p)


peptides = read_tsv(_opt(BLASTP_RAW)) if _opt(BLASTP_RAW) else read_tsv(Path(PEPTIDES))
if peptides.empty:
    peptides = read_tsv(Path(PEPTIDES))

# Immunopeptide membership = peptide length within --peptide_length (replaces the
# former 'immuno' column INTEGRATE_ENGINES used to emit).
if "length" in peptides.columns and not peptides.empty:
    _pl = PEPTIDE_LENGTH.split("-")
    peptides["immuno"] = peptides["length"].between(int(_pl[0]), int(_pl[-1]))

netmhcpan = read_tsv(_opt(NETMHCPAN_RAW))
if not netmhcpan.empty and "peptide" in netmhcpan.columns:
    peptides = peptides.merge(
        netmhcpan[["peptide", "mhc", "rank", "bind_level"]],
        on="peptide", how="left",
    )

sections = []
n = len(peptides)
n_immuno = int(peptides.get("immuno", pd.Series(dtype=bool)).astype(bool).sum())
summary = f'<p>Identified <strong>{n}</strong> peptides ' \
          f'(<strong>{n_immuno}</strong> immunopeptides).</p>'

sections.append(section("Length distribution", plot_length_histogram(peptides)))
sections.append(section("Identifications per run", plot_psms_per_run(peptides)))

if not netmhcpan.empty:
    sections.append(section("MHC class I binding (netMHCpan)",
                            plot_netmhcpan_overview(netmhcpan, "netMHCpan overview")))

mhc_ii = read_tsv(_opt(NETMHCIIPAN_RAW))
if not mhc_ii.empty:
    sections.append(section("MHC class II binding (netMHCIIpan)",
                            plot_netmhcpan_overview(mhc_ii, "netMHCIIpan overview")))

gibbs = read_tsv(_opt(GIBBS_RAW))
if not gibbs.empty:
    sections.append(section("GibbsCluster motifs", plot_gibbs_logos(gibbs)))

flashlfq = read_tsv(_opt(FLASHLFQ_RAW))
if not flashlfq.empty:
    sections.append(section("FlashLFQ quantification", quant_summary(flashlfq)))

html = HTML_TEMPLATE.format(
    sample=SAMPLE, summary=summary,
    sections="\\n".join(sections),
)
Path(OUT).write_text(html)
print(f"[report] wrote {OUT}", file=sys.stderr)

with open("versions.yml", "w") as _vf:
    _vf.write(f'"{PROCESS_NAME}":\\n')
    _vf.write(f"    python: {sys.version.split()[0]}\\n")
    _vf.write(f"    pandas: {pd.__version__}\\n")
    _vf.write(f"    matplotlib: {matplotlib.__version__}\\n")
    try:
        _vf.write(f"    logomaker: {lm.__version__}\\n")
    except (NameError, AttributeError):
        _vf.write('    logomaker: "not installed"\\n')
