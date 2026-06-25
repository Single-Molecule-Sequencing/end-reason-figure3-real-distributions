#!/usr/bin/env python3
"""Regenerate Figure 3 summary visualization from deposited Table 2 CSV."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--table",
        type=Path,
        default=Path("3_results/tables/fig3_table2_single_experiment.csv"),
        help="Deposited Table 2 CSV path",
    )
    p.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("3_results/figures/fig3_real_distributions_from_deposited_table"),
        help="Output figure prefix (without extension)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    table = args.table
    out_prefix = args.out_prefix

    if not table.exists():
        raise FileNotFoundError(f"Missing table: {table}")

    df = pd.read_csv(table)
    df = df[df["end_reason_class"] != "All focus classes"].copy()
    df["class_fraction_pct"] = pd.to_numeric(df["class_fraction_pct"], errors="coerce")
    df["mean_length_bp"] = pd.to_numeric(df["mean_length_bp"], errors="coerce")
    df["mean_qscore"] = pd.to_numeric(df["mean_qscore"], errors="coerce")

    order = ["Signal positive", "Unblock mux change", "Mux change", "Signal negative"]
    df["end_reason_class"] = pd.Categorical(df["end_reason_class"], categories=order, ordered=True)
    df = df.sort_values("end_reason_class")

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)
    colors = ["#1b9e77", "#377eb8", "#f28e2b", "#d62728"]

    axes[0].bar(df["end_reason_class"].astype(str), df["class_fraction_pct"], color=colors)
    axes[0].set_title("Fraction of focus reads (%)")
    axes[0].tick_params(axis="x", rotation=25)

    axes[1].bar(df["end_reason_class"].astype(str), df["mean_length_bp"], color=colors)
    axes[1].set_title("Mean read length (bp)")
    axes[1].tick_params(axis="x", rotation=25)

    axes[2].bar(df["end_reason_class"].astype(str), df["mean_qscore"], color=colors)
    axes[2].axhline(10, ls="--", lw=1, color="#333333")
    axes[2].set_title("Mean per-read Q-score")
    axes[2].tick_params(axis="x", rotation=25)

    run_id = str(df["run_id"].iloc[0]) if "run_id" in df.columns and not df.empty else "unknown"
    fig.suptitle(f"Figure 3 deposited-table reproducibility view ({run_id})", fontsize=12)

    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(out_prefix.with_suffix(ext), bbox_inches="tight")
    plt.close(fig)

    lineage_path = out_prefix.with_suffix(".lineage.json")
    lineage = {
        "artifact_id": out_prefix.name,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "2_analysis/scripts/fig3_from_deposited_table.py",
        "inputs": [{"path": str(table), "sha256": _sha256(table)}],
        "outputs": [
            {"path": str(out_prefix.with_suffix(".png")), "sha256": _sha256(out_prefix.with_suffix(".png"))},
            {"path": str(out_prefix.with_suffix(".pdf")), "sha256": _sha256(out_prefix.with_suffix(".pdf"))},
            {"path": str(out_prefix.with_suffix(".svg")), "sha256": _sha256(out_prefix.with_suffix(".svg"))},
        ],
        "scope": "deposited-table reproducibility figure",
    }
    lineage_path.write_text(json.dumps(lineage, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out_prefix.with_suffix('.png')}")
    print(f"wrote {out_prefix.with_suffix('.pdf')}")
    print(f"wrote {out_prefix.with_suffix('.svg')}")
    print(f"wrote {lineage_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
