#!/usr/bin/env bash
#SBATCH --job-name=sashimi_manifest
#SBATCH -n 2
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=./%x_%j.log
#SBATCH --error=./%x_%j.err

set -euo pipefail

############################################################
# Manifest-driven sashimi plotting
#
# Merged design:
#   1. Keep cohort/B1/B2/label/color loop from sashimi-20260406.sh
#   2. Use manifest-driven one-event-per-output structure
#
# Required manifest columns:
#   plot_id, cohort, AS_type, event_display, rmats_event_file
#
# Output:
#   <manifest_root>/output/<cohort>/<AS_type>/<plot_id>/
############################################################

# ── Project-level settings ────────────────────────────────────────────────
date=20260428
project_dir_name="20260428"

analysis_root="."

# Must match CONFIG$output_root in the R integrated pipeline
output_root_name="rmats_integrated_output"

manifest_root="${analysis_root}/${output_root_name}/sashimi_manifest"
manifest_tsv="${manifest_root}/sashimi_manifest_filter.tsv"

lnk_dir="./lnk"

# ── Cohort/BAM/label/color definitions from sashimi-20260406.sh ──────────
cohort_dirs=("HDAC8_MIG" "HDAC8_OE" "HDAC8_OE_vs_MIG")
cohort_labels=("HDAC8_MIG" "HDAC8_OE" "HDAC8_OE_vs_MIG")

B1s=("MIG_b1_NIR" "HDAC8OE_b1_NIR" "MIG_b1_NIR")
B2s=("MIG_b2_IR"  "HDAC8OE_b2_IR"  "HDAC8OE_b1_NIR")

L1s=("HDAC8_MIG_NIR" "HDAC8_OE_NIR" "HDAC8_MIG_NIR")
L2s=("HDAC8_MIG_IR"  "HDAC8_OE_IR"  "HDAC8_OE_NIR")

# Control/NIR color per cohort
color_Ctrl=("#00BFC4" "#FFA040" "#00BFC4")
colors_treatment=("#0B5394" "#B856D7" "#FFA040")

# Group indices from sashimi-20260406.sh
GROUP1_INDICES="1,2,3"
GROUP2_INDICES="4,5,6"

# ── Preflight checks ──────────────────────────────────────────────────────
if [[ ! -f "${manifest_tsv}" ]]; then
  echo "ERROR: manifest not found: ${manifest_tsv}" >&2
  echo "Run the R pipeline first to generate sashimi_manifest.tsv." >&2
  exit 1
fi

if ! command -v rmats2sashimiplot >/dev/null 2>&1; then
  echo "ERROR: rmats2sashimiplot not found in PATH" >&2
  exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 not found in PATH" >&2
  exit 1
fi

output_root="${manifest_root}/output"
mkdir -p "${output_root}"

# ── Main cohort loop ──────────────────────────────────────────────────────
for i in "${!cohort_labels[@]}"; do
  cohort_label="${cohort_labels[$i]}"
  cohort_dir="${cohort_dirs[$i]}"

  B1="${lnk_dir}/${B1s[$i]}.txt"
  B2="${lnk_dir}/${B2s[$i]}.txt"
  L1="${L1s[$i]}"
  L2="${L2s[$i]}"
  COLOR="${color_Ctrl[$i]},${colors_treatment[$i]}"

  if [[ ! -f "${B1}" ]]; then
    echo "WARN: B1 file missing for ${cohort_label}: ${B1}; skipping this cohort." >&2
    continue
  fi

  if [[ ! -f "${B2}" ]]; then
    echo "WARN: B2 file missing for ${cohort_label}: ${B2}; skipping this cohort." >&2
    continue
  fi

  cohort_manifest_root="${manifest_root}/by_cohort/${cohort_label}"
  mkdir -p "${cohort_manifest_root}"

  group_file="${cohort_manifest_root}/group.gf"
  {
    echo "${L1}: ${GROUP1_INDICES}"
    echo "${L2}: ${GROUP2_INDICES}"
  } > "${group_file}"

  echo "============================================================"
  echo "Cohort: ${cohort_label}"
  echo "Cohort dir: ${cohort_dir}"
  echo "B1: ${B1}"
  echo "B2: ${B2}"
  echo "Labels: ${L1} vs ${L2}"
  echo "Color: ${COLOR}"
  echo "Group file: ${group_file}"
  echo "============================================================"

  export MANIFEST_TSV="${manifest_tsv}"
  export CURRENT_COHORT="${cohort_label}"
  export B1 B2 L1 L2 COLOR
  export GROUP_FILE="${group_file}"
  export OUTPUT_ROOT="${output_root}"
  export ANALYSIS_ROOT="${analysis_root}"

  python3 - <<'PY'
import csv
import os
import subprocess
import sys

manifest_tsv = os.environ["MANIFEST_TSV"]
current_cohort = os.environ["CURRENT_COHORT"]

B1 = os.environ["B1"]
B2 = os.environ["B2"]
L1 = os.environ["L1"]
L2 = os.environ["L2"]
COLOR = os.environ["COLOR"]
group_file = os.environ["GROUP_FILE"]
output_root = os.environ["OUTPUT_ROOT"]

required_cols = ["plot_id", "cohort", "AS_type", "event_display", "rmats_event_file"]

n_total = 0
n_this_cohort = 0
n_run = 0
n_skip = 0

with open(manifest_tsv, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    missing = [c for c in required_cols if c not in reader.fieldnames]
    if missing:
        sys.exit("ERROR: manifest is missing required columns: " + ", ".join(missing))

    for row in reader:
        n_total += 1

        if row["cohort"] != current_cohort:
            continue

        n_this_cohort += 1

        plot_id = row["plot_id"]
        as_type = row["AS_type"]
        event_display = row["event_display"]
        event_file = row["rmats_event_file"]

        if not os.path.isabs(event_file):
           event_file = os.path.join(os.environ["ANALYSIS_ROOT"], event_file)

        event_file = os.path.abspath(event_file)

        if not os.path.isfile(event_file):
            print(f"WARN: event file missing for {current_cohort} {plot_id}: {event_file}; skipping", file=sys.stderr)
            n_skip += 1
            continue

        outdir = os.path.join(output_root, current_cohort, as_type, plot_id)
        os.makedirs(outdir, exist_ok=True)

        print(f"Running cohort={current_cohort} plot_id={plot_id} AS_type={as_type} event_display={event_display}")

        cmd = [
            "rmats2sashimiplot",
            "--b1", B1,
            "--b2", B2,
            "--event-type", as_type,
            "-e", event_file,
            "--group-info", group_file,
            "--font-size", "12",
            "--color", COLOR,
            "--l1", L1,
            "--l2", L2,
            "--exon_s", "1",
            "--intron_s", "5",
            "-o", outdir,
        ]

        subprocess.run(cmd, check=True)
        n_run += 1

print(
    f"Finished cohort={current_cohort}: "
    f"n_total_manifest_rows={n_total}, "
    f"n_rows_for_this_cohort={n_this_cohort}, "
    f"n_run={n_run}, "
    f"n_skip={n_skip}"
)
PY

done

echo "All cohort-specific manifest rows processed."
echo "Output root: ${output_root}"
