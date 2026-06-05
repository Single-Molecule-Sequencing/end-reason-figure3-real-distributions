#!/usr/bin/env python3
"""Figure 3 — Read-length and Q-score distributions by end-reason class.

FOR EXPERIMENTS BASECALLED WITH DORADO v1.3.1 OR LATER.
End reasons are read from the BAM er:Z: aux tag, which dorado first wrote
in v1.3.1 (DOR-1307, commit cb8fe64a, backported to release-v1.3 branch).
dorado v1.3.0 does NOT write this tag — use fig3_real_distributions_pod5_end_reasons.py
for experiments basecalled with v1.3.0 or earlier.

Both read length (query_length) and Q-score (qs tag) are also taken from
the BAM, so POD5 files are not needed.

Produces a 2×2 figure:
  top-left:     Overall read-length distribution (all reads, single black KDE)
  top-right:    Overall Q-score distribution (all reads, single black KDE)
  bottom-left:  Read length stratified by end-reason (4 colored KDEs)
  bottom-right: Q-score stratified by end-reason (4 colored KDEs)

Usage:
    python fig3_real_distributions_bam_er_tag.py \\
        --run-dir /nfs/turbo/.../20250519_1041_MN48328_AYJ384_c3faa658 \\
        --out-dir ../../3_results/figures/raw_output \\
        --peak-bp 5500   # optional: expected physical fragment peak(s) in bp
"""
from __future__ import annotations

import argparse
import hashlib
import json
import socket
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


# ---------------------------------------------------------------------------
# Display constants
# ---------------------------------------------------------------------------

FOCUS_CLASSES = [
    "signal_positive",
    "unblock_mux_change",
    "mux_change",
    "signal_negative",
]

CLASS_LABELS = {
    "signal_positive": "Signal positive",
    "unblock_mux_change": "Unblock mux change",
    "mux_change": "Mux change",
    "signal_negative": "Signal negative",
}

CLASS_COLORS = {
    "signal_positive": "#1b9e77",
    "unblock_mux_change": "#377eb8",
    "mux_change": "#f28e2b",
    "signal_negative": "#d62728",
}

LENGTH_XLIM = (1.8, 4.7)
LENGTH_TICKS = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
LENGTH_LABELS = ["100", "316", "1 kb", "3.2 kb", "10 kb", "32 kb"]
QSCORE_XLIM = (0, 25)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def read_bam_metrics(bam_dir: Path) -> pd.DataFrame:
    """Return DataFrame(read_id, read_length, qscore, end_reason) from BAM.

    end_reason is read from the er:Z: aux tag (dorado v1.3.1+).
    qscore is read from the qs tag.
    read_length is read.query_length.
    """
    bam_files = sorted(bam_dir.glob("*.bam"))
    if not bam_files:
        raise FileNotFoundError(f"No *.bam files in {bam_dir}")

    rows: list[dict] = []
    qs_missing = 0
    er_missing = 0

    for path in bam_files:
        print(f"  BAM: {path.name}", flush=True)
        with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue
                length = read.query_length
                if not length:
                    continue

                # Q-score from qs tag (dorado per-read mean quality)
                if read.has_tag("qs"):
                    qscore = float(read.get_tag("qs"))
                else:
                    qs_missing += 1
                    quals = read.query_qualities
                    if quals is None:
                        continue
                    qscore = float(np.mean(quals))

                # End reason from er:Z: tag (dorado v1.3.1+)
                if read.has_tag("er"):
                    end_reason = read.get_tag("er")
                else:
                    er_missing += 1
                    end_reason = None

                rows.append({
                    "read_id": read.query_name,
                    "read_length": length,
                    "qscore": qscore,
                    "end_reason": end_reason,
                })

    if qs_missing:
        print(f"  WARNING: {qs_missing:,} reads had no qs tag; used mean base quality",
              flush=True)
    if er_missing:
        pct = 100 * er_missing / max(len(rows), 1)
        print(f"  WARNING: {er_missing:,} reads ({pct:.1f}%) had no er:Z: tag.",
              flush=True)
        print("  Are you sure this BAM was produced by dorado v1.3.1 or later?",
              flush=True)

    data = pd.DataFrame(rows).dropna(subset=["read_length", "qscore"])
    data = data[data["read_length"] > 0].reset_index(drop=True)
    print(f"  BAM total: {len(data):,} reads", flush=True)
    return data


def build_dataset(run_dir: Path) -> pd.DataFrame:
    """Load BAM metrics (length, qscore, end_reason) for one run directory."""
    bam_dir = run_dir / "bam_pass"
    print(f"\nLoading BAM metrics from {bam_dir}", flush=True)
    data = read_bam_metrics(bam_dir)

    no_er = data["end_reason"].isna().sum()
    if no_er / len(data) > 0.05:
        raise ValueError(
            f"{no_er:,}/{len(data):,} reads missing er:Z: tag (>5%). "
            "This BAM was likely not produced by dorado v1.3.1+."
        )

    data["log10_length"] = np.log10(data["read_length"].astype(float))
    return data.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def set_plot_style() -> None:
    sns.set_theme(
        context="paper",
        style="whitegrid",
        rc={
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#333333",
            "axes.linewidth": 0.8,
            "grid.color": "#d9d9d9",
            "grid.linewidth": 0.5,
            "savefig.dpi": 300,
            "figure.dpi": 120,
        },
    )


def _shade_latest(ax: plt.Axes, color: str, alpha: float) -> None:
    if not ax.lines:
        return
    line = ax.lines[-1]
    x, y = line.get_data()
    if len(x) == 0:
        return
    ax.fill_between(x, y, 0, color=color, alpha=alpha, linewidth=0,
                    zorder=max(line.get_zorder() - 0.2, 0))
    line.set_zorder(line.get_zorder() + 0.2)


def _apply_length_axis(ax: plt.Axes) -> None:
    ax.set_xlim(*LENGTH_XLIM)
    ax.set_xticks(LENGTH_TICKS)
    ax.set_xticklabels(LENGTH_LABELS)
    ax.set_xlabel("Read length (log-scaled bp)")


def _finalize(ax: plt.Axes) -> None:
    ax.set_ylabel("KDE density")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_overall_length(ax: plt.Axes, data: pd.DataFrame, peaks: list[int]) -> None:
    if len(data) >= 5:
        sns.kdeplot(data=data, x="log10_length", ax=ax,
                    color="#333333", linewidth=1.5,
                    bw_adjust=0.85, clip=LENGTH_XLIM, warn_singular=False)
        _shade_latest(ax, "#333333", 0.14)
    for peak in peaks:
        ax.axvline(np.log10(peak), color="#222222", linestyle="--", linewidth=1.0, alpha=0.75)
        ax.text(np.log10(peak) + 0.03, 0.97, "Expected molecular length",
                transform=ax.get_xaxis_transform(), va="top", fontsize=7)
    _apply_length_axis(ax)
    _finalize(ax)
    ax.set_title("Overall read-length distribution", loc="left", fontweight="bold")
    if ax.get_legend():
        ax.get_legend().remove()


def plot_overall_qscore(ax: plt.Axes, data: pd.DataFrame) -> None:
    if len(data) >= 5:
        sns.kdeplot(data=data, x="qscore", ax=ax,
                    color="#333333", linewidth=1.5,
                    bw_adjust=0.95, clip=QSCORE_XLIM, warn_singular=False)
        _shade_latest(ax, "#333333", 0.14)
    ax.axvline(10, color="#222222", linestyle="--", linewidth=1.0, alpha=0.75)
    ax.text(10.2, 0.97, "Q10 field filter",
            transform=ax.get_xaxis_transform(), va="top", fontsize=7)
    ax.set_xlim(*QSCORE_XLIM)
    ax.set_xlabel("Mean per-read Q-score")
    _finalize(ax)
    ax.set_title("Overall Q-score distribution", loc="left", fontweight="bold")
    if ax.get_legend():
        ax.get_legend().remove()


def plot_stratified_length(ax: plt.Axes, data: pd.DataFrame, peaks: list[int]) -> None:
    for er in FOCUS_CLASSES:
        sub = data[data["end_reason"] == er]
        if len(sub) < 5:
            continue
        sns.kdeplot(data=sub, x="log10_length", ax=ax,
                    color=CLASS_COLORS[er], label=CLASS_LABELS[er],
                    linewidth=1.6 if er == "signal_positive" else 1.1,
                    bw_adjust=0.85, common_norm=False, clip=LENGTH_XLIM, warn_singular=False)
        _shade_latest(ax, CLASS_COLORS[er], 0.16 if er == "signal_positive" else 0.10)
    for peak in peaks:
        ax.axvline(np.log10(peak), color="#222222", linestyle="--", linewidth=1.0, alpha=0.75)
        ax.text(np.log10(peak) + 0.03, 0.97, "expected molecular length",
                transform=ax.get_xaxis_transform(), va="top", fontsize=7)
    _apply_length_axis(ax)
    _finalize(ax)
    ax.set_title("Read length stratified by end-reason", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="upper left")


def plot_stratified_qscore(ax: plt.Axes, data: pd.DataFrame) -> None:
    for er in FOCUS_CLASSES:
        sub = data[data["end_reason"] == er]
        if len(sub) < 5:
            continue
        sns.kdeplot(data=sub, x="qscore", ax=ax,
                    color=CLASS_COLORS[er], label=CLASS_LABELS[er],
                    linewidth=1.6 if er == "signal_positive" else 1.1,
                    bw_adjust=0.95, common_norm=False, clip=QSCORE_XLIM, warn_singular=False)
        _shade_latest(ax, CLASS_COLORS[er], 0.16 if er == "signal_positive" else 0.10)
    ax.axvline(10, color="#222222", linestyle="--", linewidth=1.0, alpha=0.75)
    ax.text(10.2, 0.97, "Q10 field filter",
            transform=ax.get_xaxis_transform(), va="top", fontsize=7)
    ax.set_xlim(*QSCORE_XLIM)
    ax.set_xlabel("Mean per-read Q-score")
    _finalize(ax)
    ax.set_title("Q-score stratified by end-reason", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="upper left")


# ---------------------------------------------------------------------------
# Main figure — 2×2 layout
# ---------------------------------------------------------------------------

def make_figure(data: pd.DataFrame, peaks: list[int], run_dir: Path, out_dir: Path) -> None:
    set_plot_style()

    focus = data[data["end_reason"].isin(FOCUS_CLASSES)].copy()

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)

    plot_overall_length(axes[0, 0], focus, peaks)
    plot_overall_qscore(axes[0, 1], focus)
    plot_stratified_length(axes[1, 0], focus, peaks)
    plot_stratified_qscore(axes[1, 1], focus)

    fig.suptitle(
        "Read-length and Q-score Distributions Stratified By End Reason",
        fontsize=12, fontweight="bold",
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "fig3_real_distributions"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {stem}.pdf", flush=True)
    print(f"Saved: {stem}.png", flush=True)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_lineage(out_dir: Path, run_dir: Path, peaks: list[int],
                  data: pd.DataFrame) -> None:
    focus = data[data["end_reason"].isin(FOCUS_CLASSES)]
    lineage = {
        "figure_id": "fig3_real_distributions",
        "paper_id": "end-reason",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "2_analysis/scripts/fig3_real_distributions_bam_er_tag.py",
        "script_sha256": _sha256(Path(__file__)),
        "host": socket.gethostname(),
        "params": {
            "run_dir": str(run_dir),
            "end_reason_source": "BAM er:Z: tag (dorado v1.3.1+, DOR-1307)",
            "qscore_source": "BAM qs tag (dorado per-read mean quality)",
            "read_length_source": "BAM query_length",
            "expected_peaks_bp": peaks,
            "qscore_reference_line": "Q10",
        },
        "read_counts_by_end_reason": {
            er: int((focus["end_reason"] == er).sum()) for er in FOCUS_CLASSES
        },
        "total_focus_reads": int(len(focus)),
    }
    (out_dir / "lineage.json").write_text(json.dumps(lineage, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--run-dir", type=Path,
        default=Path("/nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data"
                     "/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular"
                     "/20250519_1041_MN48328_AYJ384_c3faa658"),
        help="Run directory containing bam_pass/ subdirectory.",
    )
    ap.add_argument(
        "--out-dir", type=Path,
        default=Path(__file__).resolve().parents[2] / "3_results" / "figures" / "raw_output",
        help="Directory for output figures and provenance.",
    )
    ap.add_argument(
        "--peak-bp", type=int, nargs="*", default=[],
        metavar="BP",
        help="Expected physical fragment-size peak(s) in bp for dashed reference lines. "
             "Example: --peak-bp 5500",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    run_dir: Path = args.run_dir
    out_dir: Path = args.out_dir
    peaks: list[int] = sorted(set(args.peak_bp))

    if not run_dir.exists():
        print(f"ERROR: run directory not found: {run_dir}", flush=True)
        return 1

    print(f"Run directory   : {run_dir}", flush=True)
    print(f"Output directory: {out_dir}", flush=True)
    print(f"Expected peaks  : {peaks if peaks else '(none)'}", flush=True)

    data = build_dataset(run_dir)
    make_figure(data, peaks, run_dir, out_dir)
    write_lineage(out_dir, run_dir, peaks, data)
    print("Done.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
