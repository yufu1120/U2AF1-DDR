#!/usr/bin/env bash
#SBATCH --job-name=sashimi_manifest
#SBATCH -n 2
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=./%x_%j.log
#SBATCH --error=./%x_%j.err

RUNTIME_PROBE="./sashimi_${SLURM_JOB_ID:-manual}.started"

{
  echo "Script entered at $(date --iso-8601=seconds)"
  echo "Host: $(hostname)"
  echo "PID: $$"
  echo "Script: ${BASH_SOURCE[0]}"
} > "$RUNTIME_PROBE"

set -uo pipefail

# Force Python output to appear immediately in SLURM logs.
export PYTHONUNBUFFERED=1

echo "============================================================"
echo "Sashimi workflow started"
echo "Timestamp: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-not_running_under_slurm}"
echo "SLURM job name: ${SLURM_JOB_NAME:-unknown}"
echo "Submission directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Current directory: $(pwd)"
echo "Script: ${BASH_SOURCE[0]}"
echo "============================================================"

############################################################
# Cohort-aware, manifest-driven sashimi workflow
#
# 1. Read a run-specific manifest from run_sashimi_only()
# 2. Match each manifest cohort to its own B1/B2 lists, labels, and colors
# 3. Generate sashimi PDFs without aborting the whole job on one failed event
# 4. Write status/failure/missing-PDF reports
# 5. Build one combined summary PowerPoint
#
# Expected side-branch layout:
#   <analysis_root>/<output_root_name>/sashimi_runs/<sashimi_run_name>/
#       manifest/sashimi_manifest.tsv
#       plots/<cohort>/<AS_type>/<plot_id>/
#       reports/
#       logs/per_event/
#       ppt/
############################################################

# ==========================================================
# USER CONFIG: project and side-branch run
# ==========================================================
project_dir_name="20250422_IR_vs_NIR"
analysis_root="."
output_root_name="rmats_integrated_output-HDAC8_OE"
sashimi_run_name="DDR_gene_list_20260713"

# Derived manifest path:
# ${analysis_root}/${output_root_name}/sashimi_runs/${sashimi_run_name}/manifest/sashimi_manifest.tsv
lnk_dir="./lnk"

# Optional explicit manifest override. Leave empty to derive from settings above.
manifest_override=""

# PowerPoint builder. By default, use the file beside this shell script.
PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PPT_SCRIPT="${PIPELINE_DIR}/03_build_sashimi_ppt_from_manifest.py"
OVERWRITE=0
BUILD_PPT=1
DPI=200
HEARTBEAT_SECONDS=300
POLL_SECONDS=10

# ==========================================================
# USER CONFIG: cohort-specific BAM lists, labels, and colors
# Parallel arrays must have the same length and matching order.
# ==========================================================
cohort_labels=("MV4CD531" "MV4HDAC8")
B1s=("MV4CD531_b1_NIR" "MV4HDAC8_b1_NIR")
B2s=("MV4CD531_b2_IR"  "MV4HDAC8_b2_IR")
L1s=("MV4CD531_NIR" "MV4HDAC8_NIR")
L2s=("MV4CD531_IR"  "MV4HDAC8_IR")

color_group1="#85827f"
colors_group2=("#00BFC4" "#FFA040")

# Indices in the combined BAM order expected by rmats2sashimiplot group.gf.
GROUP1_INDICES="1,2,3"
GROUP2_INDICES="4,5,6"

usage() {
  cat <<EOF_USAGE
Usage: $0 [options]

Project/run options:
  --analysis-root DIR       Analysis project directory
  --output-root-name NAME   Existing integrated output folder name
  --sashimi-run-name NAME   Side-branch run name under sashimi_runs/
  --manifest FILE           Explicit manifest path; overrides derived path
  --lnk-dir DIR             Directory containing cohort B1/B2 list files

Workflow options:
  --ppt-script FILE         PowerPoint builder script
  --overwrite               Rerun events even when PDFs already exist
  --no-ppt                  Generate and validate PDFs only
  --dpi N                   PDF rendering DPI for PowerPoint (default: ${DPI})
  --heartbeat-seconds N     Main-log heartbeat interval per event (default: ${HEARTBEAT_SECONDS})
  --poll-seconds N          Child-process polling interval (default: ${POLL_SECONDS})
  -h, --help                Show this help

The cohort arrays, labels, colors, and group indices are edited in USER CONFIG.
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis-root) analysis_root="$2"; shift 2 ;;
    --output-root-name) output_root_name="$2"; shift 2 ;;
    --sashimi-run-name|--run-name) sashimi_run_name="$2"; shift 2 ;;
    --manifest) manifest_override="$2"; shift 2 ;;
    --lnk-dir) lnk_dir="$2"; shift 2 ;;
    --ppt-script) PPT_SCRIPT="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    --no-ppt) BUILD_PPT=0; shift ;;
    --dpi) DPI="$2"; shift 2 ;;
    --heartbeat-seconds) HEARTBEAT_SECONDS="$2"; shift 2 ;;
    --poll-seconds) POLL_SECONDS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

# Normalize core paths after CLI overrides.
if [[ -z "${lnk_dir}" ]]; then
  lnk_dir="${analysis_root}/../lnk"
fi
if [[ -n "${manifest_override}" ]]; then
  manifest_tsv="${manifest_override}"
  manifest_tsv="$(readlink -f "${manifest_tsv}")"
  run_dir="$(dirname "$(dirname "${manifest_tsv}")")"
else
  manifest_tsv="${analysis_root}/${output_root_name}/sashimi_runs/${sashimi_run_name}/manifest/sashimi_manifest.tsv"
  run_dir="$(dirname "$(dirname "${manifest_tsv}")")"
fi

reports_dir="${run_dir}/reports"
logs_dir="${run_dir}/logs/per_event"
ppt_dir="${run_dir}/ppt"
config_dir="${run_dir}/config/by_cohort"
mkdir -p "${reports_dir}" "${logs_dir}" "${ppt_dir}" "${config_dir}"

status_tsv="${reports_dir}/plot_status.tsv"
failed_tsv="${reports_dir}/failed_plots.tsv"
missing_tsv="${reports_dir}/missing_pdfs.tsv"
summary_txt="${reports_dir}/run_summary.txt"
cohort_config_tsv="${run_dir}/config/cohort_plot_config.tsv"

echo "[CHECKPOINT] Paths resolved"
echo "  analysis_root: ${analysis_root}"
echo "  manifest_tsv: ${manifest_tsv}"
echo "  run_dir: ${run_dir}"
echo "  lnk_dir: ${lnk_dir}"
echo "  reports_dir: ${reports_dir}"
echo "  logs_dir: ${logs_dir}"
echo "  ppt_dir: ${ppt_dir}"
echo "  heartbeat_seconds: ${HEARTBEAT_SECONDS}"
echo "  poll_seconds: ${POLL_SECONDS}"

# ==========================================================
# Preflight checks
# ==========================================================
echo "[CHECKPOINT] Starting preflight checks at $(date --iso-8601=seconds)"
if [[ ! -f "${manifest_tsv}" ]]; then
  echo "ERROR: manifest not found: ${manifest_tsv}" >&2
  echo "Run run_sashimi_only() first, or provide --manifest FILE." >&2
  exit 1
fi

if ! command -v rmats2sashimiplot >/dev/null 2>&1; then
  echo "ERROR: rmats2sashimiplot not found in PATH" >&2
  exit 1
fi
echo "[CHECKPOINT] Manifest exists"

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 not found in PATH" >&2
  exit 1
fi

echo "[CHECKPOINT] Required commands found"

n_cohorts=${#cohort_labels[@]}
for array_name in B1s B2s L1s L2s colors_group2; do
  eval "array_len=\${#${array_name}[@]}"
  if [[ "${array_len}" -ne "${n_cohorts}" ]]; then
    echo "ERROR: ${array_name} has ${array_len} entries; cohort_labels has ${n_cohorts}." >&2
    exit 1
  fi
done

# Detect duplicate cohort labels.
declare -A seen_cohorts=()
for cohort in "${cohort_labels[@]}"; do
  if [[ -n "${seen_cohorts[$cohort]:-}" ]]; then
    echo "ERROR: duplicate cohort label in configuration: ${cohort}" >&2
    exit 1
  fi
  seen_cohorts[$cohort]=1
done

# Count CSV entries, ignoring whitespace and empty entries.
count_csv_items() {
  local csv="$1"
  python3 - "$csv" <<'PY'
import sys
print(len([x for x in sys.argv[1].split(',') if x.strip()]))
PY
}

expected_group1_n="$(count_csv_items "${GROUP1_INDICES}")"
expected_group2_n="$(count_csv_items "${GROUP2_INDICES}")"

echo "[CHECKPOINT] Cohort arrays validated"

# Build one explicit cohort configuration table used by the Python runner.
printf "cohort\tB1\tB2\tL1\tL2\tcolor\tgroup_file\n" > "${cohort_config_tsv}"

for i in "${!cohort_labels[@]}"; do
  cohort="${cohort_labels[$i]}"
  echo "[CHECKPOINT] Validating cohort ${cohort} at $(date --iso-8601=seconds)"
  B1="${lnk_dir}/${B1s[$i]}.txt"
  B2="${lnk_dir}/${B2s[$i]}.txt"
  L1="${L1s[$i]}"
  L2="${L2s[$i]}"
  COLOR="${color_group1},${colors_group2[$i]}"
  cohort_dir="${config_dir}/${cohort}"
  group_file="${cohort_dir}/group.gf"
  mkdir -p "${cohort_dir}"

  for bam_list in "${B1}" "${B2}"; do
    if [[ ! -s "${bam_list}" ]]; then
      echo "ERROR: BAM-list file missing or empty for ${cohort}: ${bam_list}" >&2
      exit 1
    fi
  done

  # rMATS BAM list files are normally one comma-separated line.
  B1_n="$(python3 - "${B1}" <<'PY'
from pathlib import Path
import sys
text = Path(sys.argv[1]).read_text().strip().replace('\n', ',')
print(len([x for x in text.split(',') if x.strip()]))
PY
)"
  B2_n="$(python3 - "${B2}" <<'PY'
from pathlib import Path
import sys
text = Path(sys.argv[1]).read_text().strip().replace('\n', ',')
print(len([x for x in text.split(',') if x.strip()]))
PY
)"

  if [[ "${B1_n}" -ne "${expected_group1_n}" ]]; then
    echo "ERROR: ${cohort} B1 contains ${B1_n} BAM entries, but GROUP1_INDICES contains ${expected_group1_n} indices." >&2
    exit 1
  fi
  if [[ "${B2_n}" -ne "${expected_group2_n}" ]]; then
    echo "ERROR: ${cohort} B2 contains ${B2_n} BAM entries, but GROUP2_INDICES contains ${expected_group2_n} indices." >&2
    exit 1
  fi

  printf "%s: %s\n%s: %s\n" "${L1}" "${GROUP1_INDICES}" "${L2}" "${GROUP2_INDICES}" > "${group_file}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${cohort}" "${B1}" "${B2}" "${L1}" "${L2}" "${COLOR}" "${group_file}" \
    >> "${cohort_config_tsv}"
done

echo "[CHECKPOINT] Cohort BAM lists and group files validated"
echo "[CHECKPOINT] Cohort config table: ${cohort_config_tsv}"

# ==========================================================
# Run all configured cohorts through one manifest-aware runner
# ==========================================================
export MANIFEST_TSV="${manifest_tsv}"
export COHORT_CONFIG_TSV="${cohort_config_tsv}"
export ANALYSIS_ROOT="${analysis_root}"
export RUN_DIR="${run_dir}"
export STATUS_TSV="${status_tsv}"
export FAILED_TSV="${failed_tsv}"
export MISSING_TSV="${missing_tsv}"
export SUMMARY_TXT="${summary_txt}"
export LOGS_DIR="${logs_dir}"
export OVERWRITE
export HEARTBEAT_SECONDS
export POLL_SECONDS

echo "[CHECKPOINT] Launching manifest runner at $(date --iso-8601=seconds)"
python3 -u - <<'PY'
import csv
import os
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

manifest_path = Path(os.environ["MANIFEST_TSV"]).resolve()
cohort_config_path = Path(os.environ["COHORT_CONFIG_TSV"]).resolve()
analysis_root = Path(os.environ["ANALYSIS_ROOT"]).resolve()
run_dir = Path(os.environ["RUN_DIR"]).resolve()
status_path = Path(os.environ["STATUS_TSV"])
failed_path = Path(os.environ["FAILED_TSV"])
missing_path = Path(os.environ["MISSING_TSV"])
summary_path = Path(os.environ["SUMMARY_TXT"])
logs_dir = Path(os.environ["LOGS_DIR"])
overwrite = os.environ.get("OVERWRITE", "0") == "1"
heartbeat_seconds = max(1, int(os.environ.get("HEARTBEAT_SECONDS", "300")))
poll_seconds = max(1, int(os.environ.get("POLL_SECONDS", "10")))

required_manifest = {
    "plot_id", "cohort", "AS_type", "event_display", "rmats_event_file"
}

with cohort_config_path.open(newline="") as fh:
    configs = list(csv.DictReader(fh, delimiter="\t"))

config_by_cohort = {row["cohort"]: row for row in configs}
if len(config_by_cohort) != len(configs):
    sys.exit("ERROR: duplicated cohort labels in cohort_plot_config.tsv")

with manifest_path.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = set(reader.fieldnames or [])
    missing_cols = sorted(required_manifest - fields)
    if missing_cols:
        sys.exit("ERROR: manifest is missing required columns: " + ", ".join(missing_cols))
    manifest_rows = list(reader)

print(f"[CHECKPOINT] Manifest loaded: {len(manifest_rows)} rows from {manifest_path}", flush=True)

manifest_cohorts = sorted({row["cohort"] for row in manifest_rows})
configured_cohorts = sorted(config_by_cohort)
unconfigured = sorted(set(manifest_cohorts) - set(configured_cohorts))
if unconfigured:
    sys.exit(
        "ERROR: manifest contains cohorts without shell configuration: "
        + ", ".join(unconfigured)
    )

zero_row_cohorts = sorted(set(configured_cohorts) - set(manifest_cohorts))
if zero_row_cohorts:
    print(
        "WARN: configured cohorts with zero manifest rows: "
        + ", ".join(zero_row_cohorts),
        file=sys.stderr,
    )


def resolve_existing_file(value: str) -> Path:
    path = Path(value)
    candidates = [path] if path.is_absolute() else [analysis_root / path, run_dir / path]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return candidates[0].resolve()


def resolve_output_dir(row: dict) -> Path:
    value = (row.get("plot_output_dir") or "").strip()
    if value:
        path = Path(value)
        if path.is_absolute():
            return path
        return (analysis_root / path).resolve()
    return (run_dir / "plots" / row["cohort"] / row["AS_type"] / row["plot_id"]).resolve()


status_rows = []
total = len(manifest_rows)
for index, row in enumerate(manifest_rows, start=1):
    cohort = row["cohort"]
    cfg = config_by_cohort[cohort]
    plot_id = row["plot_id"]
    as_type = row["AS_type"]
    gene = row.get("geneSymbol", "")
    event_file = resolve_existing_file(row["rmats_event_file"])
    outdir = resolve_output_dir(row)
    outdir.mkdir(parents=True, exist_ok=True)

    safe_cohort = "".join(c if c.isalnum() or c in "._-" else "_" for c in cohort)
    safe_plot = "".join(c if c.isalnum() or c in "._-" else "_" for c in plot_id)
    log_file = logs_dir / f"{safe_cohort}__{safe_plot}.log"

    existing_pdfs = sorted(outdir.rglob("*.pdf"))
    workflow_status = "pending"
    return_code = ""
    error_message = ""

    if existing_pdfs and not overwrite:
        workflow_status = "skipped_existing"
        print(f"[{index}/{total}] SKIP existing PDFs: cohort={cohort} plot_id={plot_id}", flush=True)
    elif not event_file.is_file():
        workflow_status = "failed_missing_event_file"
        error_message = str(event_file)
        print(f"[{index}/{total}] ERROR missing event file: {event_file}", flush=True)
    else:
        cmd = [
            "rmats2sashimiplot",
            "--b1", cfg["B1"],
            "--b2", cfg["B2"],
            "--event-type", as_type,
            "-e", str(event_file),
            "--group-info", cfg["group_file"],
            "--font-size", "12",
            "--color", cfg["color"],
            "--l1", cfg["L1"],
            "--l2", cfg["L2"],
            "--exon_s", "1",
            "--intron_s", "5",
            "-o", str(outdir),
        ]

        start_time = datetime.now()
        print(
            f"[{index}/{total}] START cohort={cohort} plot_id={plot_id} "
            f"gene={gene} AS_type={as_type} "
            f"time={start_time.isoformat(timespec='seconds')} "
            f"log={log_file}",
            flush=True,
        )

        with log_file.open("w") as log:
            log.write("COMMAND:\n" + " ".join(cmd) + "\n\n")
            log.flush()

            proc = subprocess.Popen(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
            )

            last_heartbeat = time.time()
            while proc.poll() is None:
                time.sleep(poll_seconds)
                now = time.time()
                if now - last_heartbeat >= heartbeat_seconds:
                    elapsed = datetime.now() - start_time
                    try:
                        log_size = log_file.stat().st_size
                    except OSError:
                        log_size = -1
                    try:
                        pdf_count_live = len(list(outdir.rglob("*.pdf")))
                    except OSError:
                        pdf_count_live = -1
                    print(
                        f"[{index}/{total}] HEARTBEAT cohort={cohort} "
                        f"plot_id={plot_id} elapsed={elapsed} "
                        f"event_log_bytes={log_size} pdf_count={pdf_count_live}",
                        flush=True,
                    )
                    last_heartbeat = now

        end_time = datetime.now()
        elapsed = end_time - start_time
        return_code = str(proc.returncode)
        print(
            f"[{index}/{total}] END cohort={cohort} plot_id={plot_id} "
            f"return_code={proc.returncode} elapsed={elapsed}",
            flush=True,
        )

        if proc.returncode == 0:
            workflow_status = "completed"
        else:
            workflow_status = "failed_command"
            error_message = (
                f"rmats2sashimiplot exit code {proc.returncode}; see {log_file}"
            )

    pdfs = sorted(outdir.rglob("*.pdf"))
    if workflow_status in {"completed", "skipped_existing"}:
        if not pdfs:
            workflow_status = "no_pdf"
            error_message = "Command completed but no PDF was found"
        elif len(pdfs) > 1:
            workflow_status = "completed_multiple_pdfs"

    status_rows.append(
        {
            **row,
            "plot_output_dir": str(outdir),
            "resolved_rmats_event_file": str(event_file),
            "group1_label_used": cfg["L1"],
            "group2_label_used": cfg["L2"],
            "B1_used": cfg["B1"],
            "B2_used": cfg["B2"],
            "color_used": cfg["color"],
            "group_file_used": cfg["group_file"],
            "workflow_status": workflow_status,
            "n_pdfs": str(len(pdfs)),
            "pdf_files": ";".join(str(pdf) for pdf in pdfs),
            "return_code": return_code,
            "error_message": error_message,
            "log_file": str(log_file),
            "checked_at": datetime.now().isoformat(timespec="seconds"),
        }
    )

fieldnames = list(status_rows[0].keys()) if status_rows else []
subsets = {
    status_path: status_rows,
    failed_path: [
        row for row in status_rows
        if row["workflow_status"].startswith("failed")
    ],
    missing_path: [
        row for row in status_rows
        if row["workflow_status"] == "no_pdf"
    ],
}
for path, rows in subsets.items():
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

status_counts = Counter(row["workflow_status"] for row in status_rows)
cohort_counts = Counter(row["cohort"] for row in status_rows)
with summary_path.open("w") as fh:
    fh.write(f"manifest: {manifest_path}\n")
    fh.write(f"run_dir: {run_dir}\n")
    fh.write(f"total_events: {len(status_rows)}\n")
    fh.write("raw_rmats_difference: group1 minus group2\n")
    fh.write("reported_delta_psi: group2 minus group1\n")
    fh.write("\ncohort_configuration:\n")
    for cohort in configured_cohorts:
        cfg = config_by_cohort[cohort]
        fh.write(
            f"  {cohort}: group1={cfg['L1']} ({cfg['B1']}), "
            f"group2={cfg['L2']} ({cfg['B2']}), "
            f"reported_delta={cfg['L2']}_minus_{cfg['L1']}\n"
        )
    fh.write("\nmanifest_rows_by_cohort:\n")
    for cohort in sorted(cohort_counts):
        fh.write(f"  {cohort}: {cohort_counts[cohort]}\n")
    fh.write("\nworkflow_status:\n")
    for status in sorted(status_counts):
        fh.write(f"  {status}: {status_counts[status]}\n")
    if zero_row_cohorts:
        fh.write("\nconfigured_cohorts_with_zero_rows: " + ", ".join(zero_row_cohorts) + "\n")

print(f"Status table: {status_path}")
print(f"Failure table: {failed_path}")
print(f"Missing-PDF table: {missing_path}")
print(f"Summary: {summary_path}")
PY
runner_status=$?

if [[ "${runner_status}" -ne 0 ]]; then
  echo "ERROR: sashimi runner failed during manifest/config validation." >&2
  exit "${runner_status}"
fi

# ==========================================================
# Build one combined PowerPoint after all cohorts are processed
# ==========================================================
echo "[CHECKPOINT] Manifest runner completed at $(date --iso-8601=seconds)"

if [[ "${BUILD_PPT}" -eq 1 ]]; then
  echo "[CHECKPOINT] Building PowerPoint at $(date --iso-8601=seconds)"
  if [[ ! -f "${PPT_SCRIPT}" ]]; then
    echo "ERROR: PowerPoint script not found: ${PPT_SCRIPT}" >&2
    exit 1
  fi
  ppt_file="${ppt_dir}/$(basename "${run_dir}")_sashimi_summary.pptx"
  python3 -u "${PPT_SCRIPT}" \
    --manifest "${status_tsv}" \
    --out "${ppt_file}" \
    --dpi "${DPI}" \
    --skip-missing
fi

echo "[CHECKPOINT] Workflow completed at $(date --iso-8601=seconds)"
echo "Sashimi workflow finished."
echo "Manifest: ${manifest_tsv}"
echo "Run directory: ${run_dir}"
echo "Status report: ${status_tsv}"
