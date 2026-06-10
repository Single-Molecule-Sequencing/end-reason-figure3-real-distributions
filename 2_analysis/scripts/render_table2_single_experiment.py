#!/usr/bin/env python3
"""Render the source-of-record Table 2 for the Figure 3 single experiment.

The paper now uses only the Cutting-resistant E Regular run shown in Figure 3.
This script merges the two table concepts that used to live in the paper repo:

1. run-level validation metadata (POD5 reads, joined fraction, signal-positive
   prevalence, physical-window concordance, mean signal-positive length), and
2. end-reason class summary statistics (counts, fraction, length, Q-score).

Outputs:
  3_results/tables/fig3_table2_single_experiment.csv
  3_results/tables/fig3_table2_single_experiment.tex
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STATS_CSV = ROOT / "3_results" / "figures" / "fig3_summary_stats.csv"
OUT_DIR = ROOT / "3_results" / "tables"
OUT_CSV = OUT_DIR / "fig3_table2_single_experiment.csv"
OUT_TEX = OUT_DIR / "fig3_table2_single_experiment.tex"
OUT_LINEAGE = OUT_DIR / "fig3_table2_single_experiment.lineage.json"

RUN = {
    "experiment_family": "Cutting-resistant E",
    "condition": "Regular",
    "run_id": "20250519_1041_MN48328_AYJ384_c3faa658",
    "pod5_reads": "380,581",
    "joined_pct": "99.72",
    "signal_positive_pct": "96.65",
    "physical_window_pct": "88.7",
    "signal_positive_mean_bp": "4,507",
}

DISPLAY_LABEL = {
    "Signal positive": "Signal positive",
    "Unblock mux change": "Unblock mux change",
    "Mux change": "Mux change",
    "Signal negative": "Signal negative",
    "All (focus classes)": "All focus classes",
}


def _fmt_int(value: str) -> str:
    return f"{int(float(value)):,}"


def _fmt_pct(value: str) -> str:
    return f"{float(value) * 100:.2f}"


def _fmt_1(value: str) -> str:
    return f"{float(value):,.1f}"


def _tex_escape(value: str) -> str:
    return value.replace("\\", r"\textbackslash{}").replace("_", r"\_")


def read_stats() -> list[dict[str, str]]:
    with STATS_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return rows


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(rows: list[dict[str, str]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = [
        "experiment_family",
        "condition",
        "run_id",
        "pod5_reads",
        "joined_pct",
        "signal_positive_pct",
        "physical_window_pct",
        "signal_positive_mean_bp",
        "end_reason_class",
        "class_reads",
        "class_fraction_pct",
        "mean_length_bp",
        "sd_length_bp",
        "mean_qscore",
        "sd_qscore",
    ]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **RUN,
                    "end_reason_class": DISPLAY_LABEL.get(row["category"], row["category"]),
                    "class_reads": int(row["n_reads"]),
                    "class_fraction_pct": f"{float(row['fraction_of_total']) * 100:.4f}",
                    "mean_length_bp": row["mean_length_bp"],
                    "sd_length_bp": row["std_length_bp"],
                    "mean_qscore": row["mean_qscore"],
                    "sd_qscore": row["std_qscore"],
                }
            )


def render_tex(rows: list[dict[str, str]]) -> str:
    run_id_tex = _tex_escape(RUN["run_id"])
    body_rows = []
    for row in rows:
        label = DISPLAY_LABEL.get(row["category"], row["category"])
        body_rows.append(
            " & ".join(
                [
                    label,
                    _fmt_int(row["n_reads"]),
                    f"{_fmt_pct(row['fraction_of_total'])}\\%",
                    _fmt_1(row["mean_length_bp"]),
                    _fmt_1(row["std_length_bp"]),
                    f"{float(row['mean_qscore']):.2f}",
                    f"{float(row['std_qscore']):.2f}",
                ]
            )
            + r" \\"
        )
    body = "\n".join(body_rows)
    return rf"""\begin{{table}}[H]
\centering
\caption{{\textbf{{Single-experiment validation summary for Figure~\ref{{fig:quality}}.}} The table merges the Figure~3 run-level validation record with the per-end-reason read-length and Q-score summary for the Cutting-resistant E Regular run (\texttt{{{run_id_tex}}}). End reasons are read directly from POD5 acquisition files; read lengths and per-read mean Q-scores are joined from the basecaller summary/BAM by read identifier.}}
\label{{tab:fig3_single_experiment_summary}}
\scriptsize
\setlength{{\tabcolsep}}{{2.2pt}}
\renewcommand{{\arraystretch}}{{1.10}}
\begin{{tabular}}{{>{{\raggedright\arraybackslash}}p{{2.75cm}}rrrrrr}}
\toprule
\multicolumn{{7}}{{l}}{{\textbf{{Run-level validation record}}}} \\
\midrule
\multicolumn{{2}}{{l}}{{Experiment family}} & \multicolumn{{5}}{{l}}{{{RUN['experiment_family']}}} \\
\multicolumn{{2}}{{l}}{{Condition}} & \multicolumn{{5}}{{l}}{{{RUN['condition']}}} \\
\multicolumn{{2}}{{l}}{{Run ID}} & \multicolumn{{5}}{{l}}{{\texttt{{{run_id_tex}}}}} \\
\multicolumn{{2}}{{l}}{{POD5 reads}} & \multicolumn{{5}}{{r}}{{{RUN['pod5_reads']}}} \\
\multicolumn{{2}}{{l}}{{Joined to basecaller summary}} & \multicolumn{{5}}{{r}}{{{RUN['joined_pct']}\%}} \\
\multicolumn{{2}}{{l}}{{Signal-positive reads}} & \multicolumn{{5}}{{r}}{{{RUN['signal_positive_pct']}\%}} \\
\multicolumn{{2}}{{l}}{{Signal-positive reads in physical window}} & \multicolumn{{5}}{{r}}{{{RUN['physical_window_pct']}\%}} \\
\multicolumn{{2}}{{l}}{{Signal-positive mean length}} & \multicolumn{{5}}{{r}}{{{RUN['signal_positive_mean_bp']} bp}} \\
\midrule
\textbf{{End-reason class}} & \textbf{{Reads (\textit{{n}})}} & \textbf{{Fraction}} & \textbf{{Mean length (bp)}} & \textbf{{SD length (bp)}} & \textbf{{Mean Q}} & \textbf{{SD Q}} \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def main() -> int:
    rows = read_stats()
    write_csv(rows)
    OUT_TEX.write_text(render_tex(rows), encoding="utf-8")
    lineage = {
        "artifact_id": "fig3_table2_single_experiment",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "2_analysis/scripts/render_table2_single_experiment.py",
        "inputs": [
            {
                "path": "3_results/figures/fig3_summary_stats.csv",
                "sha256": _sha256(STATS_CSV),
            },
            {
                "path": "1_experiment/run_details.md",
                "role": "run metadata source",
            },
        ],
        "outputs": [
            {"path": "3_results/tables/fig3_table2_single_experiment.csv", "sha256": _sha256(OUT_CSV)},
            {"path": "3_results/tables/fig3_table2_single_experiment.tex", "sha256": _sha256(OUT_TEX)},
        ],
        "scope": "single Figure 3 experiment only",
        "run_id": RUN["run_id"],
    }
    OUT_LINEAGE.write_text(json.dumps(lineage, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_TEX}")
    print(f"wrote {OUT_LINEAGE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
