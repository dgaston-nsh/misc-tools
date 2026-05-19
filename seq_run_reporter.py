#!/usr/bin/env python3
"""
seq_run_reporter.py
-------------------
Crawls a sequencing test directory and reports FastQ, BAM, and VCF file sizes
at total, yearly, and monthly levels. Outputs TSV summaries and PDF plots.

Expected directory structure:
  <TestType>/
    <YEAR>/
      <MONTH>/
        <RunID>/
          FastQs/  (or fastq, fastqs, FASTQ, FASTQS, etc.)
          FinalBAMs/
          FinalVCFs/

Usage:
    python seq_run_reporter.py /path/to/TestType [--output-dir ./reports]
"""

import os
import re
import sys
import argparse
import csv
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# ── Non-standard dependencies ────────────────────────────────────────────────
try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: 'pandas' is not installed. See requirements below.")

try:
    import matplotlib
    matplotlib.use("Agg")  # headless
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.gridspec as gridspec
except ImportError:
    sys.exit("ERROR: 'matplotlib' is not installed. See requirements below.")

try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: 'numpy' is not installed. See requirements below.")

# ── Constants ────────────────────────────────────────────────────────────────
FASTQ_DIR_RE  = re.compile(r'^fastq[s]?$', re.IGNORECASE)
BAM_DIR_RE    = re.compile(r'^final.*bam[s]?$', re.IGNORECASE)
VCF_DIR_RE    = re.compile(r'^final.*vcf[s]?$', re.IGNORECASE)

FASTQ_EXT_RE  = re.compile(r'\.(fastq|fq)(\.gz)?$', re.IGNORECASE)
BAM_EXT_RE    = re.compile(r'\.(bam|cram)$', re.IGNORECASE)
VCF_EXT_RE    = re.compile(r'\.(vcf|vcf\.gz|bcf|gvcf|gvcf\.gz)$', re.IGNORECASE)

MONTH_NAMES = {
    1:"Jan", 2:"Feb", 3:"Mar", 4:"Apr", 5:"May",  6:"Jun",
    7:"Jul", 8:"Aug", 9:"Sep",10:"Oct",11:"Nov",  12:"Dec"
}

COLORS = {
    "fastq": "#4C72B0",
    "bam":   "#DD8452",
    "vcf":   "#55A868",
}

# ── Helpers ──────────────────────────────────────────────────────────────────

def fmt_bytes(n: float) -> str:
    """Human-readable file size."""
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(n) < 1024.0:
            return f"{n:,.2f} {unit}"
        n /= 1024.0
    return f"{n:,.2f} EB"


def classify_dir(name: str):
    """Return 'fastq', 'bam', 'vcf', or None for a subdirectory name."""
    if FASTQ_DIR_RE.match(name):
        return "fastq"
    if BAM_DIR_RE.match(name):
        return "bam"
    if VCF_DIR_RE.match(name):
        return "vcf"
    return None


def classify_file(name: str):
    """Return file category by extension, or None."""
    if FASTQ_EXT_RE.search(name):
        return "fastq"
    if BAM_EXT_RE.search(name):
        return "bam"
    if VCF_EXT_RE.search(name):
        return "vcf"
    return None


def dir_size_by_category(path: Path) -> dict:
    """
    Walk a directory and sum file sizes, classifying each file by
    extension.  Returns {category: bytes}.
    """
    totals = defaultdict(int)
    for root, _, files in os.walk(path):
        for f in files:
            cat = classify_file(f)
            if cat:
                try:
                    totals[cat] += os.path.getsize(os.path.join(root, f))
                except OSError:
                    pass
    return dict(totals)


def is_year(name: str) -> bool:
    return bool(re.fullmatch(r'\d{4}', name))


def _parse_month(name: str) -> int:
    """
    Try to extract a month number (1-12) from a directory name.
    Handles formats including (case-insensitive, any separator):
      01-JAN   02-FEB   03-MAR  ...
      JAN-01   FEB      january
      01       1        12
    Returns 0 if the name cannot be interpreted as a month.
    """
    # Split on any non-alphanumeric character so "01-JAN", "01_JAN", "01 JAN" all work
    parts = re.split(r'[^a-zA-Z0-9]+', name.strip())
    parts = [p for p in parts if p]  # drop empty strings

    month_num = 0

    for part in parts:
        # Pure number part
        if re.fullmatch(r'\d{1,2}', part):
            n = int(part)
            if 1 <= n <= 12:
                month_num = n  # keep going — a text part may confirm/override
        # Alphabetic part — try as abbreviated or full month name
        elif re.fullmatch(r'[a-zA-Z]+', part):
            try:
                m = datetime.strptime(part[:3].capitalize(), "%b").month
                return m   # text match is authoritative
            except ValueError:
                pass       # not a month name, ignore

    return month_num


def is_month(name: str) -> bool:
    """Return True if the directory name looks like a month in any supported format."""
    return _parse_month(name) != 0


def month_to_int(name: str) -> int:
    """Convert a month directory name to 1-12 (returns 0 on failure)."""
    return _parse_month(name)

# ── Core crawl ───────────────────────────────────────────────────────────────

def crawl(test_dir: Path, verbose: bool = False) -> list[dict]:
    """
    Returns a flat list of records, one per sequencing run, e.g.:
      {year, month, run, fastq_bytes, bam_bytes, vcf_bytes}
    """
    records = []

    if not test_dir.is_dir():
        sys.exit(f"ERROR: '{test_dir}' is not a directory.")

    for year_entry in sorted(test_dir.iterdir()):
        if not year_entry.is_dir() or not is_year(year_entry.name):
            continue
        year = int(year_entry.name)

        for month_entry in sorted(year_entry.iterdir()):
            if not month_entry.is_dir() or not is_month(month_entry.name):
                continue
            month = month_to_int(month_entry.name)

            for run_entry in sorted(month_entry.iterdir()):
                if not run_entry.is_dir():
                    continue
                run_name = run_entry.name
                run_totals = defaultdict(int)

                # Look inside the run directory for categorised sub-dirs AND
                # also scan loose files directly in the run dir by extension
                for sub in run_entry.iterdir():
                    if sub.is_dir():
                        cat = classify_dir(sub.name)
                        if cat:
                            # Walk the subdirectory and tally by extension
                            # (matches fastq.gz, .bam, .vcf.gz, etc.)
                            for root, _, files in os.walk(sub):
                                for f in files:
                                    fcat = classify_file(f)
                                    if fcat:
                                        try:
                                            run_totals[fcat] += os.path.getsize(
                                                os.path.join(root, f))
                                        except OSError:
                                            pass
                        else:
                            # Unknown subdir – still scan by extension
                            for root, _, files in os.walk(sub):
                                for f in files:
                                    fcat = classify_file(f)
                                    if fcat:
                                        try:
                                            run_totals[fcat] += os.path.getsize(
                                                os.path.join(root, f))
                                        except OSError:
                                            pass
                    elif sub.is_file():
                        fcat = classify_file(sub.name)
                        if fcat:
                            try:
                                run_totals[fcat] += sub.stat().st_size
                            except OSError:
                                pass

                record = {
                    "year":        year,
                    "month":       month,
                    "run":         run_name,
                    "fastq_bytes": run_totals.get("fastq", 0),
                    "bam_bytes":   run_totals.get("bam",   0),
                    "vcf_bytes":   run_totals.get("vcf",   0),
                }
                records.append(record)
                if verbose:
                    print(f"  {year}/{month:02d}/{run_name}: "
                          f"FastQ={fmt_bytes(record['fastq_bytes'])}  "
                          f"BAM={fmt_bytes(record['bam_bytes'])}  "
                          f"VCF={fmt_bytes(record['vcf_bytes'])}")

    return records

# ── Aggregation ───────────────────────────────────────────────────────────────

def build_summaries(records: list[dict], test_name: str):
    """Return (total_df, yearly_df, monthly_df) as pandas DataFrames."""
    df = pd.DataFrame(records)

    if df.empty:
        print("WARNING: No sequencing data found under the specified directory.")
        empty = pd.DataFrame()
        return empty, empty, empty

    # ── Monthly summary ──────────────────────────────────────────────────────
    monthly = (
        df.groupby(["year", "month"])[["fastq_bytes", "bam_bytes", "vcf_bytes"]]
        .sum()
        .reset_index()
    )
    monthly["total_bytes"] = monthly[["fastq_bytes", "bam_bytes", "vcf_bytes"]].sum(axis=1)
    monthly["run_count"]   = (
        df.groupby(["year", "month"])["run"].count().values
    )
    monthly["month_label"] = monthly["month"].map(MONTH_NAMES)
    monthly["period"]      = (
        monthly["year"].astype(str) + "-"
        + monthly["month"].astype(str).str.zfill(2)
    )

    # ── Yearly summary ───────────────────────────────────────────────────────
    yearly = (
        monthly.groupby("year")
        .agg(
            fastq_bytes  = ("fastq_bytes",  "sum"),
            bam_bytes    = ("bam_bytes",    "sum"),
            vcf_bytes    = ("vcf_bytes",    "sum"),
            total_bytes  = ("total_bytes",  "sum"),
            months_with_data = ("month",    "nunique"),
            run_count    = ("run_count",    "sum"),
        )
        .reset_index()
    )
    yearly["avg_fastq_per_month"] = yearly["fastq_bytes"] / yearly["months_with_data"]
    yearly["avg_bam_per_month"]   = yearly["bam_bytes"]   / yearly["months_with_data"]
    yearly["avg_vcf_per_month"]   = yearly["vcf_bytes"]   / yearly["months_with_data"]
    yearly["avg_total_per_month"] = yearly["total_bytes"] / yearly["months_with_data"]

    # ── Grand total ──────────────────────────────────────────────────────────
    total = pd.DataFrame([{
        "test_name":         test_name,
        "years_with_data":   monthly["year"].nunique(),
        "months_with_data":  len(monthly),
        "total_runs":        df["run"].count(),
        "fastq_bytes":       df["fastq_bytes"].sum(),
        "bam_bytes":         df["bam_bytes"].sum(),
        "vcf_bytes":         df["vcf_bytes"].sum(),
        "total_bytes":       df[["fastq_bytes", "bam_bytes", "vcf_bytes"]].sum().sum(),
    }])

    return total, yearly, monthly


# ── TSV output ────────────────────────────────────────────────────────────────

def write_tsvs(total_df, yearly_df, monthly_df, out_dir: Path, test_name: str):
    """Write three TSV files."""

    def add_human(df: pd.DataFrame, byte_cols: list[str]) -> pd.DataFrame:
        d = df.copy()
        for col in byte_cols:
            if col in d.columns:
                d[col + "_human"] = d[col].apply(fmt_bytes)
        return d

    byte_cols = ["fastq_bytes", "bam_bytes", "vcf_bytes", "total_bytes",
                 "avg_fastq_per_month", "avg_bam_per_month",
                 "avg_vcf_per_month",  "avg_total_per_month"]

    # Total
    t_path = out_dir / f"{test_name}_total_summary.tsv"
    add_human(total_df, byte_cols).to_csv(t_path, sep="\t", index=False)
    print(f"  Written: {t_path}")

    # Yearly
    y_path = out_dir / f"{test_name}_yearly_summary.tsv"
    add_human(yearly_df, byte_cols).to_csv(y_path, sep="\t", index=False)
    print(f"  Written: {y_path}")

    # Monthly
    m_path = out_dir / f"{test_name}_monthly_summary.tsv"
    add_human(monthly_df, byte_cols).to_csv(m_path, sep="\t", index=False)
    print(f"  Written: {m_path}")

    return t_path, y_path, m_path


# ── Plotting ──────────────────────────────────────────────────────────────────

def gb(b: float) -> float:
    return b / (1024 ** 3)


def _bar_group(ax, x, data: dict, width=0.25, labels=None):
    """Draw grouped bars for fastq/bam/vcf."""
    n = len(next(iter(data.values())))
    offsets = np.linspace(-(len(data)-1)*width/2, (len(data)-1)*width/2, len(data))
    bars = []
    for (key, vals), off in zip(data.items(), offsets):
        b = ax.bar(x + off, vals, width, label=key.upper(),
                   color=COLORS[key], alpha=0.85, edgecolor="white", linewidth=0.5)
        bars.append(b)
    if labels:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    return bars


def make_plots(total_df, yearly_df, monthly_df, out_dir: Path, test_name: str):
    pdf_path = out_dir / f"{test_name}_plots.pdf"

    with PdfPages(pdf_path) as pdf:

        # ── Page 1: Grand total pie + bar ───────────────────────────────────
        fig = plt.figure(figsize=(11, 7))
        fig.suptitle(f"{test_name}  –  Overall Storage Summary", fontsize=14, fontweight="bold")
        gs  = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)

        row = total_df.iloc[0]
        sizes  = [row["fastq_bytes"], row["bam_bytes"], row["vcf_bytes"]]
        clrs   = [COLORS["fastq"], COLORS["bam"], COLORS["vcf"]]
        labels = [f"FastQ\n{fmt_bytes(sizes[0])}", f"BAM\n{fmt_bytes(sizes[1])}", f"VCF\n{fmt_bytes(sizes[2])}"]

        ax_pie = fig.add_subplot(gs[0])
        wedges, texts, autotexts = ax_pie.pie(
            sizes, labels=labels, colors=clrs,
            autopct=lambda p: f"{p:.1f}%" if p > 0 else "",
            startangle=140, pctdistance=0.75,
            wedgeprops=dict(edgecolor="white", linewidth=1.5)
        )
        for at in autotexts:
            at.set_fontsize(9)
        ax_pie.set_title("Share by file type", fontsize=11)

        ax_bar = fig.add_subplot(gs[1])
        cats   = ["FastQ", "BAM", "VCF", "Total"]
        vals   = [gb(sizes[0]), gb(sizes[1]), gb(sizes[2]), gb(sum(sizes))]
        bars   = ax_bar.bar(cats, vals, color=[COLORS["fastq"], COLORS["bam"], COLORS["vcf"], "#9B59B6"],
                            alpha=0.85, edgecolor="white", linewidth=0.5)
        ax_bar.set_ylabel("Size (GB)", fontsize=10)
        ax_bar.set_title("Total storage by category", fontsize=11)
        ax_bar.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.1f}"))
        for bar, val in zip(bars, vals):
            ax_bar.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(vals)*0.01,
                        fmt_bytes(val * 1024**3), ha="center", va="bottom", fontsize=8)

        stats_text = (
            f"Years with data: {int(row['years_with_data'])}\n"
            f"Months with data: {int(row['months_with_data'])}\n"
            f"Total runs: {int(row['total_runs'])}\n"
            f"Grand total: {fmt_bytes(row['total_bytes'])}"
        )
        fig.text(0.5, 0.01, stats_text, ha="center", fontsize=9,
                 bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.7))

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── Page 2: Yearly grouped bar (totals) ─────────────────────────────
        if not yearly_df.empty:
            fig, axes = plt.subplots(2, 1, figsize=(11, 10))
            fig.suptitle(f"{test_name}  –  Yearly Summary", fontsize=14, fontweight="bold")

            years = yearly_df["year"].astype(str).tolist()
            x     = np.arange(len(years))

            # Top: absolute totals
            ax = axes[0]
            data = {
                "fastq": [gb(v) for v in yearly_df["fastq_bytes"]],
                "bam":   [gb(v) for v in yearly_df["bam_bytes"]],
                "vcf":   [gb(v) for v in yearly_df["vcf_bytes"]],
            }
            _bar_group(ax, x, data, labels=years)
            ax.set_ylabel("Size (GB)", fontsize=10)
            ax.set_title("Annual storage totals by file type", fontsize=11)
            ax.legend(fontsize=9)
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.1f}"))

            # Annotate months
            for i, (_, row) in enumerate(yearly_df.iterrows()):
                ax.text(i, ax.get_ylim()[1] * 0.97 if ax.get_ylim()[1] else 1,
                        f"{int(row['months_with_data'])} mo",
                        ha="center", va="top", fontsize=7, color="#555")

            # Bottom: avg per month
            ax2 = axes[1]
            data2 = {
                "fastq": [gb(v) for v in yearly_df["avg_fastq_per_month"]],
                "bam":   [gb(v) for v in yearly_df["avg_bam_per_month"]],
                "vcf":   [gb(v) for v in yearly_df["avg_vcf_per_month"]],
            }
            _bar_group(ax2, x, data2, labels=years)
            ax2.set_ylabel("Size (GB)", fontsize=10)
            ax2.set_title("Average monthly storage per year by file type", fontsize=11)
            ax2.legend(fontsize=9)
            ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.1f}"))

            fig.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ── Page 3: Monthly stacked area / line chart (all time) ────────────
        if not monthly_df.empty:
            mdf = monthly_df.sort_values(["year", "month"]).copy()
            mdf["period_idx"] = range(len(mdf))

            fig, ax = plt.subplots(figsize=(max(11, len(mdf)*0.35), 6))
            fig.suptitle(f"{test_name}  –  Monthly Storage Over Time", fontsize=14, fontweight="bold")

            x    = mdf["period_idx"].values
            fq   = [gb(v) for v in mdf["fastq_bytes"]]
            bam  = [gb(v) for v in mdf["bam_bytes"]]
            vcf  = [gb(v) for v in mdf["vcf_bytes"]]

            ax.stackplot(x, fq, bam, vcf,
                         labels=["FastQ", "BAM", "VCF"],
                         colors=[COLORS["fastq"], COLORS["bam"], COLORS["vcf"]],
                         alpha=0.80)

            # X-axis labels: show every label or thin them for readability
            tick_step = max(1, len(mdf) // 24)
            ax.set_xticks(x[::tick_step])
            ax.set_xticklabels(mdf["period"].values[::tick_step], rotation=60, ha="right", fontsize=7)
            ax.set_ylabel("Size (GB)", fontsize=10)
            ax.set_title("Stacked monthly storage (all time)", fontsize=11)
            ax.legend(loc="upper left", fontsize=9)
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.1f}"))

            # Shade year bands
            yr_boundaries = mdf.groupby("year")["period_idx"].agg(["min","max"])
            for i, (yr, row) in enumerate(yr_boundaries.iterrows()):
                ax.axvspan(row["min"] - 0.5, row["max"] + 0.5,
                           alpha=0.05 if i % 2 == 0 else 0.12,
                           color="grey", zorder=0)
                ax.text((row["min"] + row["max"]) / 2, ax.get_ylim()[1] if ax.get_ylim()[1] else 1,
                        str(yr), ha="center", va="top", fontsize=8, color="#444", alpha=0.7)

            fig.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ── Pages 4+: Per-year monthly breakdown ────────────────────────────
        if not monthly_df.empty:
            for year, ydf in monthly_df.groupby("year"):
                ydf = ydf.sort_values("month")
                fig, axes = plt.subplots(1, 2, figsize=(13, 6))
                fig.suptitle(f"{test_name}  –  Monthly Breakdown  {year}  "
                             f"({ydf['month'].nunique()} months, "
                             f"{int(ydf['run_count'].sum())} runs)",
                             fontsize=13, fontweight="bold")

                xlabels = [MONTH_NAMES[m] for m in ydf["month"]]
                x       = np.arange(len(xlabels))

                # Left: grouped bars
                ax = axes[0]
                data = {
                    "fastq": [gb(v) for v in ydf["fastq_bytes"]],
                    "bam":   [gb(v) for v in ydf["bam_bytes"]],
                    "vcf":   [gb(v) for v in ydf["vcf_bytes"]],
                }
                _bar_group(ax, x, data, labels=xlabels)
                ax.set_ylabel("Size (GB)", fontsize=10)
                ax.set_title("Monthly totals by file type", fontsize=11)
                ax.legend(fontsize=9)
                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.1f}"))

                # Right: stacked bar
                ax2 = axes[1]
                fq  = [gb(v) for v in ydf["fastq_bytes"]]
                bam = [gb(v) for v in ydf["bam_bytes"]]
                vcf = [gb(v) for v in ydf["vcf_bytes"]]
                ax2.bar(x, fq, label="FastQ", color=COLORS["fastq"], alpha=0.85, edgecolor="white")
                ax2.bar(x, bam, bottom=fq, label="BAM",
                        color=COLORS["bam"], alpha=0.85, edgecolor="white")
                ax2.bar(x, vcf, bottom=[f+b for f,b in zip(fq,bam)], label="VCF",
                        color=COLORS["vcf"], alpha=0.85, edgecolor="white")
                ax2.set_xticks(x)
                ax2.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=8)
                ax2.set_ylabel("Size (GB)", fontsize=10)
                ax2.set_title("Stacked monthly totals", fontsize=11)
                ax2.legend(fontsize=9)
                ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.1f}"))

                # Annotate run counts above each stacked bar
                for i, (_, mrow) in enumerate(ydf.iterrows()):
                    total_h = gb(mrow["fastq_bytes"] + mrow["bam_bytes"] + mrow["vcf_bytes"])
                    ax2.text(i, total_h + max(gb(ydf["fastq_bytes"]+ydf["bam_bytes"]+ydf["vcf_bytes"]))*0.01,
                             f"{int(mrow['run_count'])} runs",
                             ha="center", va="bottom", fontsize=7, color="#333")

                fig.tight_layout(rect=[0, 0, 1, 0.95])
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        # ── Metadata page ────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(8.5, 4))
        ax.axis("off")
        generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ax.text(0.5, 0.95, f"Report generated: {generated_at}",
                ha="center", va="top", fontsize=9, transform=ax.transAxes, color="#666")
        ax.text(0.5, 0.85, f"Test type: {test_name}",
                ha="center", va="top", fontsize=11, fontweight="bold", transform=ax.transAxes)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # PDF metadata
        d = pdf.infodict()
        d["Title"]   = f"{test_name} Sequencing Storage Report"
        d["Author"]  = "seq_run_reporter.py"
        d["Subject"] = "FastQ / BAM / VCF storage summary"
        d["CreationDate"] = datetime.now()

    print(f"  Written: {pdf_path}")
    return pdf_path


# ── Console summary ───────────────────────────────────────────────────────────

def print_console_summary(total_df, yearly_df, monthly_df, test_name: str):
    row = total_df.iloc[0]
    print("\n" + "═"*60)
    print(f"  TEST: {test_name}")
    print("═"*60)
    print(f"  Years with data  : {int(row['years_with_data'])}")
    print(f"  Months with data : {int(row['months_with_data'])}")
    print(f"  Total runs       : {int(row['total_runs'])}")
    print(f"  FastQ total      : {fmt_bytes(row['fastq_bytes'])}")
    print(f"  BAM total        : {fmt_bytes(row['bam_bytes'])}")
    print(f"  VCF total        : {fmt_bytes(row['vcf_bytes'])}")
    print(f"  Grand total      : {fmt_bytes(row['total_bytes'])}")
    print("═"*60)

    if not yearly_df.empty:
        print("\n  YEARLY BREAKDOWN")
        print("  " + "-"*56)
        hdr = f"  {'Year':<6} {'Months':>7} {'Runs':>6} {'FastQ':>12} {'BAM':>12} {'VCF':>12} {'Total':>12}"
        print(hdr)
        print("  " + "-"*56)
        for _, yr in yearly_df.iterrows():
            print(f"  {int(yr['year']):<6} {int(yr['months_with_data']):>7} "
                  f"{int(yr['run_count']):>6} "
                  f"{fmt_bytes(yr['fastq_bytes']):>12} "
                  f"{fmt_bytes(yr['bam_bytes']):>12} "
                  f"{fmt_bytes(yr['vcf_bytes']):>12} "
                  f"{fmt_bytes(yr['total_bytes']):>12}")
        print("\n  AVG PER ACTIVE MONTH:")
        print("  " + "-"*56)
        for _, yr in yearly_df.iterrows():
            print(f"  {int(yr['year']):<6}  FastQ={fmt_bytes(yr['avg_fastq_per_month'])}  "
                  f"BAM={fmt_bytes(yr['avg_bam_per_month'])}  "
                  f"VCF={fmt_bytes(yr['avg_vcf_per_month'])}  "
                  f"Total={fmt_bytes(yr['avg_total_per_month'])}")
    print()


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Crawl a sequencing test directory and report FastQ/BAM/VCF storage."
    )
    parser.add_argument(
        "test_dir",
        help="Path to the test directory (e.g. /data/BRCA1_Panel)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=None,
        help="Directory to write TSV and PDF reports (default: <test_dir>/reports)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print each run as it is scanned"
    )
    args = parser.parse_args()

    test_dir  = Path(args.test_dir).resolve()
    test_name = test_dir.name

    out_dir = Path(args.output_dir).resolve() if args.output_dir else test_dir / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nScanning: {test_dir}")
    print(f"Output  : {out_dir}\n")

    records = crawl(test_dir, verbose=args.verbose)

    if not records:
        print("No sequencing run data found.  Check directory structure.")
        sys.exit(1)

    total_df, yearly_df, monthly_df = build_summaries(records, test_name)

    print_console_summary(total_df, yearly_df, monthly_df, test_name)

    print("Writing TSV reports …")
    write_tsvs(total_df, yearly_df, monthly_df, out_dir, test_name)

    print("Generating plots …")
    make_plots(total_df, yearly_df, monthly_df, out_dir, test_name)

    print("\nDone.\n")


if __name__ == "__main__":
    main()
