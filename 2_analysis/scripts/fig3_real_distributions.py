#!/usr/bin/env python3
"""Figure 3 — Read-length and Q-score distributions by end-reason class.

Reads end reasons from POD5 acquisition files for a single run directory.
Joins read length and Q-score from basecalled BAM files (``qs`` tag) by read_id.
Produces a three-panel figure:
  left:   read-length KDE by end-reason (log-scaled bp)
  middle: Q-score KDE by end-reason (Q10 reference line)
  right:  % reads in expected physical fragment-size window per end-reason

Usage:
    python fig3_real_distributions.py \\
        --run-dir /nfs/turbo/.../20250519_1041_MN48328_AYJ384_c3faa658 \\
        --out-dir ../../3_results/figures/raw_output

    # optionally mark expected physical fragment peak(s) with dashed lines:
    python fig3_real_distributions.py \\
        --run-dir /nfs/turbo/.../20250519_1041_MN48328_AYJ384_c3faa658 \\
        --out-dir ../../3_results/figures/raw_output \\
        --peak-bp 3000 6000
"""
from __future__ import annotations

import argparse
import hashlib
import json
import socket
import uuid
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pod5
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
    "unblock_mux_change": "Unblock/MUX change",
    "mux_change": "MUX change",
    "signal_negative": "Signal negative",
}

CLASS_COLORS = {
    "signal_positive": "#1b9e77",
    "unblock_mux_change": "#377eb8",
    "mux_change": "#f28e2b",
    "signal_negative": "#d62728",
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _format_pod5_read_id(value: object) -> str:
    return str(uuid.UUID(bytes=bytes(value)))


def read_pod5_end_reasons(pod5_dir: Path) -> pd.DataFrame:
    """Return DataFrame(read_id, end_reason) from all *.pod5 files in pod5_dir."""
    pod5_files = sorted(pod5_dir.glob("*.pod5"))
    if not pod5_files:
        raise FileNotFoundError(f"No *.pod5 files found in {pod5_dir}")

    frames: list[pd.DataFrame] = []
    for pod5_path in pod5_files:
        print(f"  reading POD5 end reasons: {pod5_path.name}", flush=True)
        reader = pod5.Reader(str(pod5_path))
        try:
            for batch_idx in range(reader.read_table.num_record_batches):
                batch = reader.read_table.get_batch(batch_idx).select(["read_id", "end_reason"])
                chunk = batch.to_pandas()
                if chunk.empty:
                    continue
                chunk["read_id"] = [_format_pod5_read_id(v) for v in chunk["read_id"]]
                chunk["end_reason"] = chunk["end_reason"].astype(str)
                frames.append(chunk)
        finally:
            reader.close()

    if not frames:
        raise ValueError(f"No reads extracted from POD5 files in {pod5_dir}")

    data = pd.concat(frames, ignore_index=True)
    # resolve duplicate read_ids (keep first occurrence)
    dupes = data.groupby("read_id")["end_reason"].nunique()
    conflicts = dupes[dupes > 1]
    if not conflicts.empty:
        raise ValueError(f"Conflicting end reasons for {len(conflicts)} read IDs in {pod5_dir}")
    data = data.drop_duplicates(subset=["read_id"], keep="first")
    print(f"  POD5: {len(data):,} reads with end reasons", flush=True)
    return data[["read_id", "end_reason"]].reset_index(drop=True)


def read_bam_metrics(bam_dir: Path) -> pd.DataFrame:
    """Return DataFrame(read_id, read_length, qscore) from all BAM files.

    Q-score is taken from the ``qs`` tag (dorado's per-read mean quality score).
    Falls back to mean base quality from ``query_qualities`` if ``qs`` is absent.
    Read length is ``query_length`` (full query sequence including soft clips).
    """
    bam_files = sorted(bam_dir.glob("*.bam"))
    if not bam_files:
        raise FileNotFoundError(f"No *.bam files found in {bam_dir}")

    rows: list[dict] = []
    qs_missing = 0
    for bam_path in bam_files:
        print(f"  reading BAM metrics: {bam_path.name}", flush=True)
        with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue
                read_id = read.query_name
                read_length = read.query_length
                if read_length is None or read_length == 0:
                    continue
                # prefer qs tag; fall back to mean base quality
                if read.has_tag("qs"):
                    qscore = float(read.get_tag("qs"))
                else:
                    qs_missing += 1
                    quals = read.query_qualities
                    if quals is None:
                        continue
                    qscore = float(np.mean(quals))
                rows.append({"read_id": read_id, "read_length": read_length, "qscore": qscore})

    if not rows:
        raise ValueError(f"No reads extracted from BAM files in {bam_dir}")
    if qs_missing:
        print(f"  WARNING: {qs_missing:,} reads lacked qs tag; used mean base quality instead", flush=True)

    data = pd.DataFrame(rows)
    data = data.dropna(subset=["read_length", "qscore"])
    data = data[data["read_length"] > 0]
    print(f"  BAM: {len(data):,} reads with length + Q-score", flush=True)
    return data.reset_index(drop=True)


def build_dataset(run_dir: Path) -> pd.DataFrame:
    """Join POD5 end reasons with BAM metrics for a single run directory."""
    pod5_dir = run_dir / "pod5"
    bam_dir = run_dir / "bam_pass"

    print(f"Loading POD5 end reasons from {pod5_dir}", flush=True)
    er = read_pod5_end_reasons(pod5_dir)

    print(f"Loading BAM metrics (qs tag) from {bam_dir}", flush=True)
    bam = read_bam_metrics(bam_dir)

    merged = bam.merge(er, on="read_id", how="inner")
    n_unmatched = len(bam) - len(merged)
    if len(bam) > 0 and n_unmatched / len(bam) > 0.05:
        raise ValueError(
            f"More than 5% of BAM reads ({n_unmatched:,}/{len(bam):,}) "
            "had no matching POD5 end reason"
        )
    print(f"  Joined: {len(merged):,} reads ({n_unmatched:,} BAM reads unmatched in POD5)", flush=True)

    merged = merged[merged["end_reason"].isin(FOCUS_CLASSES)].copy()
    merged["log10_length"] = np.log10(merged["read_length"].astype(float))
    merged["end_reason_label"] = merged["end_reason"].map(CLASS_LABELS)
    return merged.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Physical-window helpers
# ---------------------------------------------------------------------------

def in_physical_window(lengths: np.ndarray, peaks: list[int]) -> np.ndarray:
    if not peaks:
        return np.zeros(len(lengths), dtype=bool)
    mask = np.zeros(len(lengths), dtype=bool)
    for center in peaks:
        tol = max(50.0, 0.075 * float(center))
        mask |= np.abs(lengths - center) <= tol
    return mask


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


def plot_length_kde(ax: plt.Axes, data: pd.DataFrame, peaks: list[int], show_legend: bool) -> None:
    xlim = (1.8, 4.7)
    for er in FOCUS_CLASSES:
        sub = data[data["end_reason"] == er]
        if len(sub) < 5:
            continue
        sns.kdeplot(data=sub, x="log10_length", ax=ax,
                    color=CLASS_COLORS[er], label=CLASS_LABELS[er],
                    linewidth=1.6 if er == "signal_positive" else 1.1,
                    bw_adjust=0.85, common_norm=False, clip=xlim, warn_singular=False)
        _shade_latest(ax, CLASS_COLORS[er], 0.16 if er == "signal_positive" else 0.10)

    for peak in sorted(set(peaks)):
        ax.axvline(np.log10(peak), color="#555555", linestyle="--", linewidth=0.8, alpha=0.45)
    if peaks:
        ax.text(0.01, 0.97, "Dashed: expected physical peaks",
                transform=ax.transAxes, va="top", fontsize=7, color="#555555")

    ticks = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    labels = ["100", "316", "1 kb", "3.2 kb", "10 kb", "32 kb"]
    ax.set_xlim(*xlim)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Read length (log-scaled bp)")
    ax.set_ylabel("KDE density")
    ax.set_title("Read length by end-reason", loc="left", fontsize=9, fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_legend:
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    elif ax.get_legend():
        ax.get_legend().remove()


def plot_qscore_kde(ax: plt.Axes, data: pd.DataFrame) -> None:
    xlim = (0, 26)
    for er in FOCUS_CLASSES:
        sub = data[data["end_reason"] == er]
        if len(sub) < 5:
            continue
        sns.kdeplot(data=sub, x="qscore", ax=ax,
                    color=CLASS_COLORS[er], label=CLASS_LABELS[er],
                    linewidth=1.6 if er == "signal_positive" else 1.1,
                    bw_adjust=0.95, common_norm=False, clip=xlim, warn_singular=False)
        _shade_latest(ax, CLASS_COLORS[er], 0.16 if er == "signal_positive" else 0.10)

    ax.axvline(10, color="#222222", linestyle="--", linewidth=1.0, alpha=0.75)
    ax.text(10.2, 0.96, "Q10 field filter",
            transform=ax.get_xaxis_transform(), va="top", fontsize=7)
    ax.set_xlim(*xlim)
    ax.set_xlabel("Mean per-read Q-score (qs tag)")
    ax.set_ylabel("KDE density")
    ax.set_title("Q-score by end-reason", loc="left", fontsize=9, fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ax.get_legend():
        ax.get_legend().remove()


def plot_physical_bar(ax: plt.Axes, data: pd.DataFrame, peaks: list[int]) -> None:
    rows = []
    lengths = data["read_length"].to_numpy(dtype=np.float64)
    in_window = in_physical_window(lengths, peaks)
    for er in FOCUS_CLASSES:
        mask = data["end_reason"].to_numpy() == er
        n = int(mask.sum())
        if n == 0:
            continue
        pct = 100.0 * int((in_window & mask).sum()) / n
        rows.append({"end_reason": er, "label": CLASS_LABELS[er], "pct": pct, "n": n})

    if not rows:
        ax.set_visible(False)
        return

    df = pd.DataFrame(rows)
    colors = [CLASS_COLORS[v] for v in df["end_reason"]]
    ax.barh(df["label"], df["pct"], color=colors, alpha=0.88)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Reads in expected physical window (%)")
    ax.set_title("Physical-window match", loc="left", fontsize=9, fontweight="bold")
    for idx, row in enumerate(df.itertuples(index=False)):
        ax.text(min(row.pct + 1.5, 97), idx, f"{row.pct:.1f}%", va="center", fontsize=7)
    ax.grid(axis="x", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if not peaks:
        ax.text(0.5, 0.5, "No expected peaks provided\n(all bars will be 0%)",
                transform=ax.transAxes, ha="center", va="center", fontsize=8, color="#888")


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure(data: pd.DataFrame, peaks: list[int], run_dir: Path, out_dir: Path) -> None:
    set_plot_style()
    fig, axes = plt.subplots(
        nrows=1, ncols=3,
        figsize=(13.5, 4.8),
        gridspec_kw={"width_ratios": [1.2, 1.0, 0.8]},
        constrained_layout=True,
    )
    plot_length_kde(axes[0], data, peaks, show_legend=True)
    plot_qscore_kde(axes[1], data)
    plot_physical_bar(axes[2], data, peaks)

    run_name = run_dir.name
    fig.suptitle(
        f"Signal-positive reads recapitulate known physical molecule sizes\n"
        f"Run: {run_name}",
        fontsize=11, fontweight="bold",
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fig3_real_distributions.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "fig3_real_distributions.png", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to {out_dir / 'fig3_real_distributions.pdf'}", flush=True)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_lineage(
    out_dir: Path,
    run_dir: Path,
    peaks: list[int],
    data: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    counts = {er: int((data["end_reason"] == er).sum()) for er in FOCUS_CLASSES}
    lineage = {
        "figure_id": "fig3_real_distributions",
        "paper_id": "end-reason",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "2_analysis/scripts/fig3_real_distributions.py",
        "script_sha256": _sha256(Path(__file__)),
        "host": socket.gethostname(),
        "params": {
            "run_dir": str(run_dir),
            "end_reason_source": "POD5 read_table end_reason field",
            "qscore_source": "BAM qs tag (dorado per-read mean quality)",
            "read_length_source": "BAM query_length",
            "expected_peaks_bp": peaks,
            "qscore_reference_line": "Q10",
        },
        "read_counts_by_end_reason": counts,
        "total_focus_reads": int(len(data)),
    }
    (out_dir / "lineage.json").write_text(json.dumps(lineage, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--run-dir", type=Path,
        default=Path("/nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data"
                     "/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular"
                     "/20250519_1041_MN48328_AYJ384_c3faa658"),
        help="Run directory containing pod5/ and bam_pass/ subdirectories.",
    )
    ap.add_argument(
        "--out-dir", type=Path,
        default=Path(__file__).resolve().parents[2] / "3_results" / "figures" / "raw_output",
        help="Directory to write output figures and provenance.",
    )
    ap.add_argument(
        "--peak-bp", type=int, nargs="*", default=[],
        metavar="BP",
        help="Expected physical fragment-size peak(s) in bp. "
             "Used to draw dashed reference lines and compute physical-window percentages. "
             "Example: --peak-bp 3000 6000",
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

    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Run directory : {run_dir}", flush=True)
    print(f"Output directory: {out_dir}", flush=True)
    print(f"Expected peaks  : {peaks} bp", flush=True)

    data = build_dataset(run_dir)

    print("Plotting figure...", flush=True)
    make_figure(data, peaks, run_dir, out_dir)

    write_lineage(out_dir, run_dir, peaks, data, args)
    print("Done.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
