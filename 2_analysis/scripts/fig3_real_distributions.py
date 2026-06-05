#!/usr/bin/env python3
"""Build real-data Figure 3 distributions for the end-reason paper.

The script starts from the basecalled MinKNOW run folders for the three
proof-of-principle cohorts on Great Lakes/SMB. For each run it records the basecalled BAM/FASTQ directories, extracts
end_reason directly from POD5 acquisition records, and joins read length and
Q-score from basecaller-generated sequencing_summary_*.txt tables.

Default outputs:
  Y:\paper-data\end-reason\fig3_real_distributions\
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import html
import json
import math
import re
import socket
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pod5
import seaborn as sns


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

COHORT_LABELS = {
    "single_species_eco53ki": "Single-species Eco53KI",
    "rc_eco53ki_pvuii": "Eco53KI + PvuII",
    "cutting_res_e": "Cutting-resistant E",
}

RUN_LABELS = {
    "2925_0_5uL_25fmol_Srf1": "SrfI 2925",
    "2925_100fmol_2uL_no2ndRE": "No second enzyme",
    "2925_PvuII_0_5uL_25fmol": "PvuII 2925",
    "pCYP_2987": "pCYP 2987",
    "pCYP_3040": "pCYP 3040",
    "30_40_PVU2": "PvuII 30/40",
    "30_40_SRF1": "SrfI 30/40",
    "Regular": "Regular",
}

COHORT_NOTES = {
    "single_species_eco53ki": (
        "Defined single-molecule species controls from the Lab Updates deck; "
        "restriction fragments create known physical species before sequencing."
    ),
    "rc_eco53ki_pvuii": (
        "Eco53KI+PvuII restriction controls with two short, physically separable "
        "species in the 0.47/0.63 kb range."
    ),
    "cutting_res_e": (
        "Restriction-enzyme proof set from the Lab Updates deck; SrfI/PvuII "
        "produce separable known-source fragments and the regular prep retains "
        "the full-length plasmid species."
    ),
}

SUMMARY_USECOLS = [
    "read_id",
    "sequence_length_template",
    "mean_qscore_template",
    "end_reason",
]

BASECALLED_DIR_NAMES = ["bam_pass", "fastq_pass"]
POD5_END_REASON_TABLE = "fig3_pod5_end_reasons.parquet"


@dataclass(frozen=True)
class Cohort:
    cohort_id: str
    source_dir: Path


def default_share() -> Path:
    candidates = [
        Path(r"Y:\gregfar\SMS\SMS_POP_data"),
        Path("/nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def default_out_dir() -> Path:
    candidates = [
        Path(r"Y:\paper-data\end-reason\fig3_real_distributions"),
        Path("/nfs/turbo/umms-atheylab/paper-data/end-reason/fig3_real_distributions"),
    ]
    for candidate in candidates:
        if candidate.parent.exists():
            return candidate
    return Path(__file__).resolve().parents[1] / "data" / "processed" / "fig3_real_distributions"


def default_peak_stats() -> Path:
    candidates = [
        Path(r"Y:\paper-data\end-reason\fig3_peak_stats\per_run_peak_er.csv"),
        Path("/nfs/turbo/umms-atheylab/paper-data/end-reason/fig3_peak_stats/per_run_peak_er.csv"),
        Path(__file__).resolve().parents[1] / "data" / "processed" / "fig3_peak_stats" / "per_run_peak_er.csv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def cohorts_from_share(share: Path) -> list[Cohort]:
    return [
        Cohort("single_species_eco53ki", share / "Single_MOlecule_Species_Eco53KI"),
        Cohort("rc_eco53ki_pvuii", share / "SMS_RC_Eco53KI_PvuII"),
        Cohort("cutting_res_e", share / "Single_Molecule_Seqeuncing_Cutting_Res_E"),
    ]


def find_summaries(cohorts: Iterable[Cohort]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cohort in cohorts:
        if not cohort.source_dir.exists():
            raise FileNotFoundError(f"Missing cohort directory: {cohort.source_dir}")
        for path in sorted(cohort.source_dir.glob("*/*/sequencing_summary_*.txt")):
            run_dir = path.parent
            basecalled_dirs = {name: run_dir / name for name in BASECALLED_DIR_NAMES}
            pod5_dir = run_dir / "pod5"
            pod5_files = sorted(pod5_dir.glob("*.pod5")) if pod5_dir.exists() else []
            if not pod5_files:
                raise FileNotFoundError(f"No POD5 files found for {cohort.cohort_id}/{path.parent.parent.name}: {pod5_dir}")
            rows.append(
                {
                    "cohort": cohort.cohort_id,
                    "run_label": path.parent.parent.name,
                    "path": path,
                    "bytes": path.stat().st_size,
                    "run_dir": run_dir,
                    "bam_pass_dir": basecalled_dirs["bam_pass"],
                    "fastq_pass_dir": basecalled_dirs["fastq_pass"],
                    "bam_pass_entries": count_dir_entries(basecalled_dirs["bam_pass"]),
                    "fastq_pass_entries": count_dir_entries(basecalled_dirs["fastq_pass"]),
                    "pod5_dir": pod5_dir,
                    "pod5_files": len(pod5_files),
                    "pod5_bytes": sum(file.stat().st_size for file in pod5_files),
                }
            )
    if not rows:
        raise FileNotFoundError("No sequencing_summary_*.txt files found for Figure 3 cohorts")
    return rows


def count_dir_entries(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return sum(1 for _ in path.iterdir())
    except OSError:
        return -1


def format_pod5_read_id(value: object) -> str:
    return str(uuid.UUID(bytes=bytes(value)))


def read_pod5_end_reason_file(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    reader = pod5.Reader(str(path))
    frames: list[pd.DataFrame] = []
    counts: Counter[str] = Counter()
    try:
        for batch_index in range(reader.read_table.num_record_batches):
            batch = reader.read_table.get_batch(batch_index).select(["read_id", "end_reason"])
            chunk = batch.to_pandas()
            if chunk.empty:
                continue
            chunk["read_id"] = [format_pod5_read_id(value) for value in chunk["read_id"]]
            chunk["end_reason"] = chunk["end_reason"].astype(str)
            counts.update(chunk["end_reason"].tolist())
            frames.append(chunk)
    finally:
        reader.close()

    if frames:
        data = pd.concat(frames, ignore_index=True)
    else:
        data = pd.DataFrame(columns=["read_id", "end_reason"])
    file_row = {
        "pod5_file": str(path),
        "bytes": path.stat().st_size,
        "n_reads": int(len(data)),
        "end_reason_counts_json": json.dumps(dict(sorted(counts.items())), sort_keys=True),
    }
    return data, file_row


def extract_pod5_end_reasons(
    run_inputs: list[dict[str, object]],
    out_dir: Path,
) -> tuple[dict[tuple[str, str], pd.DataFrame], pd.DataFrame]:
    by_run: dict[tuple[str, str], pd.DataFrame] = {}
    file_rows: list[dict[str, object]] = []
    all_frames: list[pd.DataFrame] = []

    for item in run_inputs:
        cohort = str(item["cohort"])
        run_label = str(item["run_label"])
        pod5_dir = Path(item["pod5_dir"])
        run_frames: list[pd.DataFrame] = []
        print(f"reading POD5 end reasons {cohort}/{run_label}: {pod5_dir}", flush=True)
        for pod5_path in sorted(pod5_dir.glob("*.pod5")):
            pod5_data, file_row = read_pod5_end_reason_file(pod5_path)
            file_row.update({"cohort": cohort, "run_label": run_label})
            file_rows.append(file_row)
            if not pod5_data.empty:
                run_frames.append(pod5_data)

        if not run_frames:
            raise ValueError(f"No POD5 reads extracted for {cohort}/{run_label}: {pod5_dir}")

        run_data = pd.concat(run_frames, ignore_index=True)
        duplicate_end_reasons = run_data.groupby("read_id", sort=False)["end_reason"].nunique()
        conflicts = duplicate_end_reasons[duplicate_end_reasons > 1]
        if not conflicts.empty:
            raise ValueError(
                f"Conflicting POD5 end reasons for {len(conflicts)} duplicate read IDs in {cohort}/{run_label}"
            )
        run_data = run_data.drop_duplicates(subset=["read_id"], keep="first")
        run_data["end_reason"] = run_data["end_reason"].astype("category")
        by_run[(cohort, run_label)] = run_data[["read_id", "end_reason"]].copy()

        export = run_data.copy()
        export.insert(0, "run_label", run_label)
        export.insert(0, "cohort", cohort)
        all_frames.append(export)

    pod5_file_manifest = pd.DataFrame(file_rows)
    pod5_file_manifest.to_csv(out_dir / "fig3_pod5_file_manifest.csv", index=False)

    pod5_end_reasons = pd.concat(all_frames, ignore_index=True)
    pod5_end_reasons.to_parquet(out_dir / POD5_END_REASON_TABLE, index=False)
    return by_run, pod5_file_manifest


def stable_seed(*parts: object) -> int:
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()
    return int(digest[:8], 16)


def safe_slug(value: object) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "-", str(value)).strip("-").lower()
    return slug or "run"


def display_run_label(value: object) -> str:
    return RUN_LABELS.get(str(value), str(value).replace("_", " "))


def load_expected_peaks(path: Path) -> dict[str, list[int]]:
    if not path.exists():
        return {}
    peak_stats = pd.read_csv(path)
    peak_stats = peak_stats[
        (peak_stats["end_reason"] == "signal_positive")
        & (peak_stats["n_reads"] >= 1_000)
        & peak_stats["peak_center_bp"].notna()
    ]
    peaks: dict[str, list[int]] = {}
    for run_label, sub in peak_stats.groupby("run_label", sort=True):
        centers = sorted({int(v) for v in sub["peak_center_bp"].tolist() if v > 0})
        if centers:
            peaks[str(run_label)] = centers
    return peaks


def in_physical_window(lengths: np.ndarray, centers: list[int]) -> np.ndarray:
    if not centers:
        return np.zeros(lengths.size, dtype=bool)
    matched = np.zeros(lengths.size, dtype=bool)
    for center in centers:
        tolerance = max(50.0, 0.075 * float(center))
        matched |= np.abs(lengths - center) <= tolerance
    return matched


def update_sample(
    reservoirs: dict[tuple[str, str, str], pd.DataFrame],
    key: tuple[str, str, str],
    incoming: pd.DataFrame,
    target: int,
    rng: np.random.Generator,
) -> None:
    if incoming.empty:
        return
    sampled = incoming.copy()
    sampled["_sample_rank"] = rng.random(len(sampled))
    current = reservoirs.get(key)
    combined = sampled if current is None else pd.concat([current, sampled], ignore_index=True)
    if len(combined) > target:
        combined = combined.nsmallest(target, "_sample_rank")
    reservoirs[key] = combined


def update_summary_rows(
    summary_acc: dict[tuple[str, str, str], dict[str, float]],
    cohort: str,
    run_label: str,
    sub: pd.DataFrame,
    expected_peaks: dict[str, list[int]],
) -> None:
    lengths = sub["sequence_length_template"].to_numpy(dtype=np.float64, copy=False)
    qscores = sub["mean_qscore_template"].to_numpy(dtype=np.float64, copy=False)
    physical = in_physical_window(lengths, expected_peaks.get(run_label, []))
    q10 = qscores >= 10.0
    for end_reason in FOCUS_CLASSES:
        mask = sub["end_reason"].to_numpy() == end_reason
        if not mask.any():
            continue
        key = (cohort, run_label, end_reason)
        acc = summary_acc.setdefault(
            key,
            {
                "n_reads": 0,
                "q10_pass_reads": 0,
                "physical_window_reads": 0,
                "length_sum": 0.0,
            },
        )
        n = int(mask.sum())
        acc["n_reads"] += n
        acc["q10_pass_reads"] += int(q10[mask].sum())
        acc["physical_window_reads"] += int(physical[mask].sum())
        acc["length_sum"] += float(lengths[mask].sum())


def stream_summaries(
    summaries: list[dict[str, object]],
    pod5_end_reasons: dict[tuple[str, str], pd.DataFrame],
    expected_peaks: dict[str, list[int]],
    sample_per_run_class: int,
    chunksize: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reservoirs: dict[tuple[str, str, str], pd.DataFrame] = {}
    summary_acc: dict[tuple[str, str, str], dict[str, float]] = {}
    concordance: Counter[tuple[str, str, str, str]] = Counter()
    join_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(20260521)

    for item in summaries:
        cohort = str(item["cohort"])
        run_label = str(item["run_label"])
        path = Path(item["path"])
        pod5_map = pod5_end_reasons[(cohort, run_label)].rename(columns={"end_reason": "pod5_end_reason"})
        pod5_read_ids = set(pod5_map["read_id"].tolist())
        seen_summary_ids: set[str] = set()
        summary_rows = 0
        valid_summary_rows = 0
        joined_rows = 0
        missing_pod5_rows = 0
        mismatched_rows = 0
        print(f"joining summary length/Q-score to POD5 end reasons {cohort}/{run_label}: {path}", flush=True)
        for chunk in pd.read_csv(
            path,
            sep="\t",
            usecols=SUMMARY_USECOLS,
            chunksize=chunksize,
            engine="c",
        ):
            summary_rows += int(len(chunk))
            seen_summary_ids.update(chunk["read_id"].astype(str).tolist())
            chunk = chunk.rename(
                columns={
                    "sequence_length_template": "read_length",
                    "mean_qscore_template": "qscore",
                    "end_reason": "summary_end_reason",
                }
            )
            chunk = chunk[
                chunk["read_length"].notna()
                & chunk["qscore"].notna()
                & (chunk["read_length"] > 0)
                & np.isfinite(chunk["qscore"])
            ].copy()
            valid_summary_rows += int(len(chunk))
            if chunk.empty:
                continue

            chunk["read_id"] = chunk["read_id"].astype(str)
            chunk = chunk.merge(pod5_map, on="read_id", how="left")
            missing_pod5 = chunk["pod5_end_reason"].isna()
            missing_pod5_rows += int(missing_pod5.sum())
            chunk = chunk[~missing_pod5].copy()
            if chunk.empty:
                continue

            chunk["end_reason"] = chunk["pod5_end_reason"].astype(str)
            joined_rows += int(len(chunk))
            mismatches = chunk["summary_end_reason"].astype(str) != chunk["end_reason"]
            mismatched_rows += int(mismatches.sum())
            pairs = chunk.groupby(["end_reason", "summary_end_reason"], sort=False).size()
            for (pod5_reason, summary_reason), count in pairs.items():
                concordance[(cohort, run_label, str(pod5_reason), str(summary_reason))] += int(count)

            chunk = chunk[chunk["end_reason"].isin(FOCUS_CLASSES)].copy()
            if chunk.empty:
                continue

            stats_chunk = chunk.rename(
                columns={
                    "read_length": "sequence_length_template",
                    "qscore": "mean_qscore_template",
                }
            )
            update_summary_rows(summary_acc, cohort, run_label, stats_chunk, expected_peaks)

            chunk["cohort"] = cohort
            chunk["run_label"] = run_label
            chunk["log10_length"] = np.log10(chunk["read_length"].astype(float))
            chunk["physical_window"] = in_physical_window(
                chunk["read_length"].to_numpy(dtype=np.float64),
                expected_peaks.get(run_label, []),
            )
            chunk["end_reason_label"] = chunk["end_reason"].map(CLASS_LABELS)
            keep = [
                "cohort",
                "run_label",
                "end_reason",
                "end_reason_label",
                "read_length",
                "log10_length",
                "qscore",
                "physical_window",
            ]
            for end_reason, sub in chunk.groupby("end_reason", sort=False):
                key = (cohort, run_label, str(end_reason))
                update_sample(reservoirs, key, sub[keep], sample_per_run_class, rng)

        if summary_rows and missing_pod5_rows / summary_rows > 0.01:
            raise ValueError(
                f"More than 1% of summary rows lack POD5 end reasons for {cohort}/{run_label}: "
                f"{missing_pod5_rows:,}/{summary_rows:,}"
            )
        reads_in_both = len(pod5_read_ids.intersection(seen_summary_ids))
        join_rows.append(
            {
                "cohort": cohort,
                "run_label": run_label,
                "pod5_reads": len(pod5_read_ids),
                "summary_rows": summary_rows,
                "valid_summary_rows": valid_summary_rows,
                "reads_in_both": reads_in_both,
                "reads_only_in_pod5": len(pod5_read_ids) - reads_in_both,
                "reads_only_in_summary": len(seen_summary_ids) - reads_in_both,
                "joined_valid_rows": joined_rows,
                "summary_rows_missing_pod5": missing_pod5_rows,
                "summary_pod5_end_reason_mismatches": mismatched_rows,
                "summary_pod5_end_reason_mismatch_pct": 100.0 * mismatched_rows / joined_rows
                if joined_rows
                else math.nan,
            }
        )

    sample = pd.concat(reservoirs.values(), ignore_index=True)
    sample = sample.drop(columns=["_sample_rank"], errors="ignore")
    sample["cohort_label"] = sample["cohort"].map(COHORT_LABELS)

    summary_rows: list[dict[str, object]] = []
    for (cohort, run_label, end_reason), acc in sorted(summary_acc.items()):
        n_reads = int(acc["n_reads"])
        summary_rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABELS.get(cohort, cohort),
                "run_label": run_label,
                "end_reason": end_reason,
                "end_reason_label": CLASS_LABELS.get(end_reason, end_reason),
                "n_reads": n_reads,
                "q10_pass_reads": int(acc["q10_pass_reads"]),
                "q10_pass_pct": 100.0 * acc["q10_pass_reads"] / n_reads if n_reads else math.nan,
                "physical_window_reads": int(acc["physical_window_reads"]),
                "physical_window_pct": (
                    100.0 * acc["physical_window_reads"] / n_reads if n_reads else math.nan
                ),
                "mean_length": acc["length_sum"] / n_reads if n_reads else math.nan,
            }
        )
    summary = pd.DataFrame(summary_rows)
    concordance_rows = [
        {
            "cohort": cohort,
            "run_label": run_label,
            "pod5_end_reason": pod5_reason,
            "summary_end_reason": summary_reason,
            "n_reads": count,
        }
        for (cohort, run_label, pod5_reason, summary_reason), count in sorted(concordance.items())
    ]
    return sample, summary, pd.DataFrame(concordance_rows), pd.DataFrame(join_rows)


def write_input_manifests(out_dir: Path, summaries: list[dict[str, object]]) -> None:
    rows = []
    for item in summaries:
        rows.append(
            {
                "cohort": item["cohort"],
                "run_label": item["run_label"],
                "run_dir": str(item["run_dir"]),
                "sequencing_summary": str(item["path"]),
                "sequencing_summary_bytes": item["bytes"],
                "pod5_dir": str(item["pod5_dir"]),
                "pod5_files": item["pod5_files"],
                "pod5_bytes": item["pod5_bytes"],
                "bam_pass_dir": str(item["bam_pass_dir"]),
                "bam_pass_entries": item["bam_pass_entries"],
                "fastq_pass_dir": str(item["fastq_pass_dir"]),
                "fastq_pass_entries": item["fastq_pass_entries"],
                "end_reason_source": "pod5",
                "length_qscore_source": "sequencing_summary",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(out_dir / "fig3_pod5_input_manifest.csv", index=False)
    manifest.to_csv(out_dir / "fig3_basecalled_input_manifest.csv", index=False)


def write_eight_run_summary(
    out_dir: Path,
    summary: pd.DataFrame,
    summaries: list[dict[str, object]],
    join_stats: pd.DataFrame,
) -> pd.DataFrame:
    """Write one compact row per Great Lakes proof-of-principle run."""
    join_lookup = {
        (str(row.cohort), str(row.run_label)): row
        for row in join_stats.itertuples(index=False)
    } if not join_stats.empty else {}
    rows: list[dict[str, object]] = []
    for item in sorted(summaries, key=lambda x: (str(x["cohort"]), str(x["run_label"]))):
        cohort = str(item["cohort"])
        run_label = str(item["run_label"])
        run_summary = summary[(summary["cohort"] == cohort) & (summary["run_label"] == run_label)]
        signal = run_summary[run_summary["end_reason"] == "signal_positive"]
        unblock = run_summary[run_summary["end_reason"] == "unblock_mux_change"]
        signal_row = signal.iloc[0] if not signal.empty else None
        unblock_row = unblock.iloc[0] if not unblock.empty else None
        focus_reads = int(run_summary["n_reads"].sum())
        q10_pass = int(run_summary["q10_pass_reads"].sum()) if "q10_pass_reads" in run_summary else 0
        joined = join_lookup.get((cohort, run_label))
        pod5_reads = int(getattr(joined, "pod5_reads", 0)) if joined is not None else math.nan
        reads_in_both = int(getattr(joined, "reads_in_both", 0)) if joined is not None else math.nan
        rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABELS.get(cohort, cohort),
                "run_label": run_label,
                "focus_reads": focus_reads,
                "q10_pass": q10_pass,
                "physical_window_reads": (
                    int(signal_row["physical_window_reads"]) if signal_row is not None else math.nan
                ),
                "signal_positive_reads": int(signal_row["n_reads"]) if signal_row is not None else math.nan,
                "signal_positive_q10_pct": float(signal_row["q10_pass_pct"]) if signal_row is not None else math.nan,
                "signal_positive_physical_pct": (
                    float(signal_row["physical_window_pct"]) if signal_row is not None else math.nan
                ),
                "signal_positive_mean_length": float(signal_row["mean_length"]) if signal_row is not None else math.nan,
                "unblock_reads": int(unblock_row["n_reads"]) if unblock_row is not None else math.nan,
                "unblock_q10_pct": float(unblock_row["q10_pass_pct"]) if unblock_row is not None else math.nan,
                "unblock_mean_length": float(unblock_row["mean_length"]) if unblock_row is not None else math.nan,
                "pod5_reads": pod5_reads,
                "summary_rows": int(getattr(joined, "summary_rows", 0)) if joined is not None else math.nan,
                "reads_in_both": reads_in_both,
                "reads_only_in_pod5": (
                    int(getattr(joined, "reads_only_in_pod5", 0)) if joined is not None else math.nan
                ),
                "reads_only_in_summary": (
                    int(getattr(joined, "reads_only_in_summary", 0)) if joined is not None else math.nan
                ),
                "summary_pod5_end_reason_mismatch_pct": (
                    float(getattr(joined, "summary_pod5_end_reason_mismatch_pct", math.nan))
                    if joined is not None
                    else math.nan
                ),
                "pod5_files": int(item["pod5_files"]),
                "pod5_bytes": int(item["pod5_bytes"]),
                "sequencing_summary_bytes": int(item["bytes"]),
                "bam_pass_entries": int(item["bam_pass_entries"]),
                "fastq_pass_entries": int(item["fastq_pass_entries"]),
                "signal_positive_pct": (
                    100 * int(signal_row["n_reads"]) / focus_reads if signal_row is not None and focus_reads else math.nan
                ),
                "unblock_pct": (
                    100 * int(unblock_row["n_reads"]) / focus_reads if unblock_row is not None and focus_reads else math.nan
                ),
                "joined_pct_of_pod5": (
                    100 * reads_in_both / pod5_reads
                    if isinstance(pod5_reads, int) and pod5_reads
                    else math.nan
                ),
            }
        )
    run_summary = pd.DataFrame(rows)
    run_summary.to_csv(out_dir / "fig3_eight_run_summary.csv", index=False)
    return run_summary


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


def shade_latest_density(ax: plt.Axes, color: str, alpha: float) -> None:
    """Fill the area under the most recently drawn KDE curve."""
    if not ax.lines:
        return
    line = ax.lines[-1]
    x_data, y_data = line.get_data()
    if len(x_data) == 0 or len(y_data) == 0:
        return
    ax.fill_between(
        x_data,
        y_data,
        0,
        color=color,
        alpha=alpha,
        linewidth=0,
        zorder=max(line.get_zorder() - 0.2, 0),
    )
    line.set_zorder(line.get_zorder() + 0.2)


def apply_length_axis(ax: plt.Axes) -> None:
    ticks = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    labels = ["100", "316", "1 kb", "3.2 kb", "10 kb", "32 kb"]
    ax.set_xlim(1.8, 4.7)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Read length (log-scaled bp)")


def plot_kde(
    ax: plt.Axes,
    data: pd.DataFrame,
    x: str,
    xlim: tuple[float, float],
    xlabel: str,
    expected_peaks: list[int] | None = None,
    show_legend: bool = False,
) -> None:
    for end_reason in FOCUS_CLASSES:
        sub = data[data["end_reason"] == end_reason]
        if len(sub) < 5:
            continue
        sns.kdeplot(
            data=sub,
            x=x,
            ax=ax,
            color=CLASS_COLORS[end_reason],
            label=CLASS_LABELS[end_reason],
            linewidth=1.6 if end_reason == "signal_positive" else 1.1,
            bw_adjust=0.85 if x == "log10_length" else 0.95,
            common_norm=False,
            clip=xlim,
            warn_singular=False,
        )
        shade_latest_density(ax, CLASS_COLORS[end_reason], 0.16 if end_reason == "signal_positive" else 0.10)
    if expected_peaks:
        for peak in sorted(set(expected_peaks)):
            ax.axvline(np.log10(peak), color="#555555", linestyle="--", linewidth=0.8, alpha=0.45)
    if x == "qscore":
        ax.axvline(10, color="#222222", linestyle="--", linewidth=1.0, alpha=0.75)
        ax.text(10.2, 0.96, "Q10 field filter", transform=ax.get_xaxis_transform(), va="top", fontsize=7)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("KDE density")
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_legend:
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    elif ax.get_legend() is not None:
        ax.get_legend().remove()


def plot_overall_kde(
    ax: plt.Axes,
    data: pd.DataFrame,
    x: str,
    xlim: tuple[float, float],
    xlabel: str,
    expected_peaks: list[int] | None = None,
) -> None:
    if len(data) >= 5:
        sns.kdeplot(
            data=data,
            x=x,
            ax=ax,
            color="#333333",
            label="All focus-class reads",
            linewidth=1.5,
            bw_adjust=0.85 if x == "log10_length" else 0.95,
            clip=xlim,
            warn_singular=False,
        )
        shade_latest_density(ax, "#333333", 0.14)
    if expected_peaks:
        for peak in sorted(set(expected_peaks)):
            ax.axvline(np.log10(peak), color="#555555", linestyle="--", linewidth=0.8, alpha=0.45)
    if x == "qscore":
        ax.axvline(10, color="#222222", linestyle="--", linewidth=1.0, alpha=0.75)
        ax.text(10.2, 0.96, "Q10 field filter", transform=ax.get_xaxis_transform(), va="top", fontsize=7)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("KDE density")
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ax.get_legend() is not None:
        ax.get_legend().remove()


def cohort_expected_peaks(expected_peaks: dict[str, list[int]], run_labels: Iterable[str]) -> list[int]:
    peaks: list[int] = []
    for run_label in run_labels:
        peaks.extend(expected_peaks.get(str(run_label), []))
    return sorted(set(peaks))


def plot_physical_bar(ax: plt.Axes, summary: pd.DataFrame, cohort: str, run_label: str | None = None) -> None:
    sub = summary[summary["cohort"] == cohort].copy()
    if run_label is not None:
        sub = sub[sub["run_label"] == run_label].copy()
    agg = (
        sub.groupby(["end_reason", "end_reason_label"], as_index=False)
        .agg(n_reads=("n_reads", "sum"), physical_window_reads=("physical_window_reads", "sum"))
    )
    agg["physical_window_pct"] = 100 * agg["physical_window_reads"] / agg["n_reads"]
    agg["order"] = agg["end_reason"].map({v: i for i, v in enumerate(FOCUS_CLASSES)})
    agg = agg.sort_values("order")
    colors = [CLASS_COLORS[v] for v in agg["end_reason"]]
    ax.barh(agg["end_reason_label"], agg["physical_window_pct"], color=colors, alpha=0.88)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Reads in expected physical windows (%)")
    ax.set_ylabel("")
    for idx, row in enumerate(agg.itertuples(index=False)):
        ax.text(
            min(float(row.physical_window_pct) + 1.5, 98),
            idx,
            f"{row.physical_window_pct:.1f}%",
            va="center",
            fontsize=7,
        )
    ax.grid(axis="x", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_combined_figure(
    sample: pd.DataFrame,
    summary: pd.DataFrame,
    expected_peaks: dict[str, list[int]],
    out_dir: Path,
) -> None:
    set_plot_style()
    fig, axes = plt.subplots(
        nrows=3,
        ncols=3,
        figsize=(12.5, 9.2),
        gridspec_kw={"width_ratios": [1.18, 1.0, 0.82]},
        constrained_layout=True,
    )
    for row, cohort in enumerate(COHORT_LABELS):
        cohort_sample = sample[sample["cohort"] == cohort]
        peaks = cohort_expected_peaks(expected_peaks, cohort_sample["run_label"].unique())
        label = COHORT_LABELS[cohort]

        ax_len = axes[row, 0]
        plot_kde(
            ax_len,
            cohort_sample,
            x="log10_length",
            xlim=(1.8, 4.7),
            xlabel="",
            expected_peaks=peaks,
            show_legend=(row == 0),
        )
        apply_length_axis(ax_len)
        ax_len.set_title(f"{label}: read length by end-reason", loc="left", fontsize=9, fontweight="bold")
        ax_len.text(
            0.01,
            0.96,
            "Dashed lines: expected physical peaks",
            transform=ax_len.transAxes,
            va="top",
            fontsize=7,
            color="#555555",
        )

        ax_q = axes[row, 1]
        plot_kde(
            ax_q,
            cohort_sample,
            x="qscore",
            xlim=(0, 26),
            xlabel="Mean per-read Q-score",
            expected_peaks=None,
            show_legend=False,
        )
        ax_q.set_title(f"{label}: Q-score by end-reason", loc="left", fontsize=9, fontweight="bold")

        ax_bar = axes[row, 2]
        plot_physical_bar(ax_bar, summary, cohort)
        ax_bar.set_title(f"{label}: physical match", loc="left", fontsize=9, fontweight="bold")

    fig.suptitle(
        "Signal-positive reads recapitulate known physical molecule sizes",
        fontsize=13,
        fontweight="bold",
    )
    fig.savefig(out_dir / "fig3_real_distributions.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "fig3_real_distributions.png", bbox_inches="tight")
    plt.close(fig)


def make_per_cohort_figures(
    sample: pd.DataFrame,
    summary: pd.DataFrame,
    expected_peaks: dict[str, list[int]],
    out_dir: Path,
) -> None:
    set_plot_style()
    for cohort, label in COHORT_LABELS.items():
        cohort_sample = sample[sample["cohort"] == cohort]
        peaks = cohort_expected_peaks(expected_peaks, cohort_sample["run_label"].unique())
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)

        plot_overall_kde(
            axes[0, 0],
            cohort_sample,
            x="log10_length",
            xlim=(1.8, 4.7),
            xlabel="",
            expected_peaks=peaks,
        )
        apply_length_axis(axes[0, 0])
        axes[0, 0].set_title("Overall read-length distribution", loc="left", fontweight="bold")

        plot_overall_kde(axes[0, 1], cohort_sample, x="qscore", xlim=(0, 26), xlabel="Mean per-read Q-score")
        axes[0, 1].set_title("Overall Q-score distribution", loc="left", fontweight="bold")

        plot_kde(
            axes[1, 0],
            cohort_sample,
            x="log10_length",
            xlim=(1.8, 4.7),
            xlabel="",
            expected_peaks=peaks,
            show_legend=True,
        )
        apply_length_axis(axes[1, 0])
        axes[1, 0].set_title("Read length stratified by end-reason", loc="left", fontweight="bold")

        plot_kde(
            axes[1, 1],
            cohort_sample,
            x="qscore",
            xlim=(0, 26),
            xlabel="Mean per-read Q-score",
            show_legend=True,
        )
        axes[1, 1].set_title("Q-score stratified by end-reason", loc="left", fontweight="bold")

        fig.suptitle(f"{label}: Figure 3 proof-of-principle distributions", fontsize=12, fontweight="bold")
        safe = cohort.replace("_", "-")
        fig.savefig(out_dir / f"fig3_{safe}_distributions.pdf", bbox_inches="tight")
        fig.savefig(out_dir / f"fig3_{safe}_distributions.png", bbox_inches="tight")
        plt.close(fig)


def make_per_run_figures(
    sample: pd.DataFrame,
    summary: pd.DataFrame,
    expected_peaks: dict[str, list[int]],
    out_dir: Path,
) -> list[dict[str, str]]:
    set_plot_style()
    run_rows: list[dict[str, str]] = []
    for (cohort, run_label), run_sample in sample.groupby(["cohort", "run_label"], sort=True):
        label = COHORT_LABELS.get(str(cohort), str(cohort))
        peaks = expected_peaks.get(str(run_label), [])
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)

        plot_overall_kde(
            axes[0, 0],
            run_sample,
            x="log10_length",
            xlim=(1.8, 4.7),
            xlabel="",
            expected_peaks=peaks,
        )
        apply_length_axis(axes[0, 0])
        axes[0, 0].set_title("Overall read-length distribution", loc="left", fontweight="bold")

        plot_overall_kde(axes[0, 1], run_sample, x="qscore", xlim=(0, 26), xlabel="Mean per-read Q-score")
        axes[0, 1].set_title("Overall Q-score distribution", loc="left", fontweight="bold")

        plot_kde(
            axes[1, 0],
            run_sample,
            x="log10_length",
            xlim=(1.8, 4.7),
            xlabel="",
            expected_peaks=peaks,
            show_legend=True,
        )
        apply_length_axis(axes[1, 0])
        axes[1, 0].set_title("Read length stratified by end-reason", loc="left", fontweight="bold")

        plot_kde(
            axes[1, 1],
            run_sample,
            x="qscore",
            xlim=(0, 26),
            xlabel="Mean per-read Q-score",
            show_legend=True,
        )
        axes[1, 1].set_title("Q-score stratified by end-reason", loc="left", fontweight="bold")

        fig.suptitle(
            f"{label} / {display_run_label(run_label)}: read-length and Q-score distributions",
            fontsize=12,
            fontweight="bold",
        )
        safe = f"{safe_slug(cohort)}_{safe_slug(run_label)}"
        pdf = f"fig3_run_{safe}_distributions.pdf"
        png = f"fig3_run_{safe}_distributions.png"
        fig.savefig(out_dir / pdf, bbox_inches="tight")
        fig.savefig(out_dir / png, bbox_inches="tight")
        plt.close(fig)
        run_rows.append(
            {
                "cohort": str(cohort),
                "cohort_label": label,
                "run_label": str(run_label),
                "pdf": pdf,
                "png": png,
            }
        )
    return run_rows


def write_viewer(out_dir: Path, sample: pd.DataFrame, summary: pd.DataFrame, summaries: list[dict[str, object]]) -> None:
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    cohort_cards = []
    for cohort, label in COHORT_LABELS.items():
        safe = cohort.replace("_", "-")
        sub = summary[summary["cohort"] == cohort]
        total_reads = int(sub["n_reads"].sum())
        signal = sub[sub["end_reason"] == "signal_positive"]
        signal_match = (
            100 * signal["physical_window_reads"].sum() / signal["n_reads"].sum()
            if not signal.empty and signal["n_reads"].sum()
            else math.nan
        )
        cohort_cards.append(
            f"""
            <section class="card view hidden" id="{safe_slug(cohort)}">
              <h2>{html.escape(label)}</h2>
              <p>{html.escape(COHORT_NOTES[cohort])}</p>
              <p><b>Total focus-class reads:</b> {total_reads:,} &nbsp; <b>Signal-positive physical-window match:</b> {signal_match:.1f}%</p>
              <img src="fig3_{safe}_distributions.png" alt="{html.escape(label)} distributions">
            </section>
            """
        )
    run_cards = []
    run_buttons = []
    for (cohort, run_label), _run_sample in sample.groupby(["cohort", "run_label"], sort=True):
        run_id = f"run-{safe_slug(cohort)}-{safe_slug(run_label)}"
        safe = f"{safe_slug(cohort)}_{safe_slug(run_label)}"
        png = f"fig3_run_{safe}_distributions.png"
        label = COHORT_LABELS.get(str(cohort), str(cohort))
        sub = summary[(summary["cohort"] == cohort) & (summary["run_label"] == run_label)]
        total_reads = int(sub["n_reads"].sum())
        signal = sub[sub["end_reason"] == "signal_positive"]
        signal_match = (
            100 * signal["physical_window_reads"].sum() / signal["n_reads"].sum()
            if not signal.empty and signal["n_reads"].sum()
            else math.nan
        )
        display_label = display_run_label(run_label)
        run_buttons.append(f'<button data-target="{run_id}">{html.escape(display_label)}</button>')
        run_cards.append(
            f"""
            <section class="card view hidden" id="{run_id}">
              <h2>{html.escape(label)} / {html.escape(display_label)}</h2>
              <p><b>Total focus-class reads:</b> {total_reads:,} &nbsp; <b>Signal-positive physical-window match:</b> {signal_match:.1f}%</p>
              <img src="{png}" alt="{html.escape(str(run_label))} per-run distributions">
            </section>
            """
        )
    sources = "\n".join(
        (
            f"<li><b>{html.escape(str(item['run_label']))}</b>: "
            f"POD5 <code>{html.escape(str(item['pod5_dir']))}</code> "
            f"({int(item['pod5_files']):,} files; {int(item['pod5_bytes']):,} bytes); "
            f"length/Q-score join <code>{html.escape(str(item['path']))}</code></li>"
        )
        for item in summaries
    )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Figure 3 POD5-derived end-reason distributions</title>
  <style>
    body {{ font-family: Arial, Helvetica, sans-serif; margin: 28px; color: #222; background: #f7f7f4; }}
    h1 {{ margin-bottom: 0.2rem; }}
    .meta {{ color: #555; margin-top: 0; }}
    .card {{ background: white; border: 1px solid #ddd; padding: 18px; margin: 18px 0; box-shadow: 0 1px 4px rgba(0,0,0,0.05); }}
    .nav {{ position: sticky; top: 0; background: #f7f7f4; padding: 8px 0; border-bottom: 1px solid #ddd; z-index: 2; }}
    button {{ border: 1px solid #aaa; background: white; padding: 7px 10px; margin: 0 6px 6px 0; cursor: pointer; }}
    button.active {{ background: #222; color: white; }}
    .hidden {{ display: none; }}
    img {{ width: 100%; max-width: 1250px; border: 1px solid #ddd; background: white; }}
    code {{ font-size: 0.9em; }}
    li {{ margin: 0.3rem 0; }}
  </style>
</head>
<body>
  <h1>Figure 3 POD5-derived end-reason distributions</h1>
  <p class="meta">Generated {html.escape(generated)} on {html.escape(socket.gethostname())}. End reasons are read directly from POD5 acquisition records. Read lengths and Q-scores are joined by read ID from basecaller-generated <code>sequencing_summary_*.txt</code> files. Q-score panels include the Q10 field-filter reference line.</p>
  <div class="nav">
    <button class="active" data-target="combined">Combined</button>
    {''.join(f'<button data-target="{safe_slug(cohort)}">{html.escape(label)}</button>' for cohort, label in COHORT_LABELS.items())}
    {''.join(run_buttons)}
    <button data-target="source">Source + provenance</button>
  </div>
  <section class="card">
    <h2>Original prompt/problem</h2>
    <p>Use the three Great Lakes proof-of-principle datasets to make Figure 3 of the end-reason paper, registering the datasets and tracking bioinformatic workflows/figure-generation scripts. Q-score graphs should include a Q10 reference line. The purpose is to show that only signal-positive reads recapitulate the known physical molecules. The current rebuild starts from the POD5 directories and uses POD5-derived <code>end_reason</code> values, not sequencing-summary-derived end reasons.</p>
    <p><b>Copilot session:</b> <a href="file:///C:/Users/gregfar/.copilot/session-state/71a931ec-a7fa-472d-99db-41b15c8d6eb9/">file:///C:/Users/gregfar/.copilot/session-state/71a931ec-a7fa-472d-99db-41b15c8d6eb9/</a></p>
  </section>
  <section class="card view" id="combined">
    <h2>Combined Figure 3 candidate</h2>
    <p>Read-length KDEs are plotted on log-scaled base-pair coordinates. Dashed length lines mark signal-positive physical peak centers from the existing trough-bounded peak-stat workflow. The right column quantifies exact read counts inside those physical windows.</p>
    <img src="fig3_real_distributions.png" alt="Combined Figure 3 candidate">
  </section>
  {''.join(cohort_cards)}
  {''.join(run_cards)}
  <section class="card view hidden" id="source">
    <h2>Source evidence and provenance</h2>
    <p>Experimental expectations were cross-checked against <code>C:\\Users\\gregfar\\Downloads\\Lab Updates.pptx</code>; extracted slide text is saved in the Copilot session files. POD5 input discovery is tracked in <code>fig3_pod5_input_manifest.csv</code> and per-file POD5 read counts are tracked in <code>fig3_pod5_file_manifest.csv</code>. End-reason provenance is the POD5 <code>end_reason</code> field; <code>sequencing_summary_*.txt</code> is used only for read length and Q-score by read-ID join. Concordance with summary end reasons is in <code>fig3_pod5_summary_concordance.csv</code>.</p>
    <ul>{sources}</ul>
  </section>
  <script>
    const buttons = document.querySelectorAll("button[data-target]");
    const views = document.querySelectorAll(".view");
    for (const button of buttons) {{
      button.addEventListener("click", () => {{
        for (const other of buttons) other.classList.toggle("active", other === button);
        for (const view of views) view.classList.toggle("hidden", view.id !== button.dataset.target);
      }});
    }}
  </script>
</body>
</html>
"""
    (out_dir / "fig3_real_distribution_viewer.html").write_text(html_text, encoding="utf-8")


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def output_manifest(out_dir: Path) -> list[dict[str, object]]:
    rows = []
    for path in sorted(out_dir.iterdir()):
        if path.is_file() and path.name != "lineage.json":
            rows.append(
                {
                    "name": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    return rows


def write_lineage(
    out_dir: Path,
    summaries: list[dict[str, object]],
    expected_peaks: dict[str, list[int]],
    sample: pd.DataFrame,
    summary: pd.DataFrame,
    join_stats: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    lineage = {
        "figure_id": "fig3_real_distributions",
        "paper_id": "end-reason",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "figures_src/fig3_real_distributions.py",
        "script_sha256": sha256_file(Path(__file__)),
        "host": socket.gethostname(),
        "params": {
            "sample_per_run_class": args.sample_per_run_class,
            "chunksize": args.chunksize,
            "reuse_existing": args.reuse_existing,
            "source_type": "pod5_end_reason_joined_to_basecalled_summary",
            "end_reason_source": "POD5 read_table end_reason field, read directly from POD5 acquisition files",
            "length_qscore_source": "sequencing_summary_*.txt joined by read_id",
            "qscore_reference_line": "Q10",
            "physical_window": "+/- max(50 bp, 7.5% of signal_positive peak center)",
        },
        "inputs": [
            {
                "cohort": item["cohort"],
                "run_label": item["run_label"],
                "run_dir": str(item["run_dir"]),
                "pod5_dir": str(item["pod5_dir"]),
                "pod5_files": item["pod5_files"],
                "pod5_bytes": item["pod5_bytes"],
                "sequencing_summary": str(item["path"]),
                "sequencing_summary_bytes": item["bytes"],
                "bam_pass_dir": str(item["bam_pass_dir"]),
                "bam_pass_entries": item["bam_pass_entries"],
                "fastq_pass_dir": str(item["fastq_pass_dir"]),
                "fastq_pass_entries": item["fastq_pass_entries"],
            }
            for item in summaries
        ],
        "expected_peaks_by_run": expected_peaks,
        "outputs": {
            "sample_rows": int(len(sample)),
            "summary_rows": int(len(summary)),
            "join_stats": join_stats.to_dict(orient="records"),
            "files": output_manifest(out_dir),
        },
        "deck_evidence": {
            "path": r"C:\Users\gregfar\Downloads\Lab Updates.pptx",
            "notes": [
                "Slide text describes defined single-molecule species preparation, restriction-enzyme cuts, and separable fragment lengths.",
                "Slide 10 states read length shows two peaks for SrfI-cut DNA and that the source sequence is known before sequencing.",
            ],
        },
    }
    (out_dir / "lineage.json").write_text(json.dumps(lineage, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--share", type=Path, default=default_share())
    ap.add_argument("--out-dir", type=Path, default=default_out_dir())
    ap.add_argument("--peak-stats", type=Path, default=default_peak_stats())
    ap.add_argument("--sample-per-run-class", type=int, default=25_000)
    ap.add_argument("--chunksize", type=int, default=750_000)
    ap.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Reuse fig3_distribution_sample.parquet and summary CSV in --out-dir instead of re-reading raw summaries.",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    cohorts = cohorts_from_share(args.share)
    summaries = find_summaries(cohorts)
    expected_peaks = load_expected_peaks(args.peak_stats)
    write_input_manifests(out_dir, summaries)
    if args.reuse_existing:
        sample = pd.read_parquet(out_dir / "fig3_distribution_sample.parquet")
        summary = pd.read_csv(out_dir / "fig3_distribution_summary.csv")
        join_stats_path = out_dir / "fig3_pod5_summary_join_stats.csv"
        join_stats = pd.read_csv(join_stats_path) if join_stats_path.exists() else pd.DataFrame()
    else:
        pod5_end_reasons, _pod5_file_manifest = extract_pod5_end_reasons(summaries, out_dir)
        sample, summary, concordance, join_stats = stream_summaries(
            summaries=summaries,
            pod5_end_reasons=pod5_end_reasons,
            expected_peaks=expected_peaks,
            sample_per_run_class=args.sample_per_run_class,
            chunksize=args.chunksize,
        )
        sample.to_parquet(out_dir / "fig3_distribution_sample.parquet", index=False)
        sample.to_csv(out_dir / "fig3_distribution_sample.csv", index=False)
        summary.to_csv(out_dir / "fig3_distribution_summary.csv", index=False)
        concordance.to_csv(out_dir / "fig3_pod5_summary_concordance.csv", index=False)
        join_stats.to_csv(out_dir / "fig3_pod5_summary_join_stats.csv", index=False)

    write_eight_run_summary(out_dir, summary, summaries, join_stats)
    make_combined_figure(sample, summary, expected_peaks, out_dir)
    make_per_cohort_figures(sample, summary, expected_peaks, out_dir)
    make_per_run_figures(sample, summary, expected_peaks, out_dir)
    write_viewer(out_dir, sample, summary, summaries)
    write_lineage(out_dir, summaries, expected_peaks, sample, summary, join_stats, args)

    print(f"wrote Figure 3 real-data outputs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
