# RNA-seq splicing analysis

Each cohort has its own self-contained scripts.

## Standard contents

- `01-rMATS/`: SLURM scripts that run rMATS for individual cohort comparisons.
- `02-downstream-analysis/01_rmats_integrated_pipeline.R`: imports, summarize, and annotates rMATS JC/JCEC results, applies count/PSI/P-value/FDR filters, makes delta PSI graphs, writes master tables, and performs splice-site sequence analysis.
- `02-downstream-analysis/02_run_sashimi_from_manifest.sh`: runs manifest-selected events through `rmats2sashimiplot`.
    - the `rmats2sashimiplot` needs to be install in the environment first.
- `02-downstream-analysis/03_build_sashimi_ppt_from_manifest.py`: assembles sashimi PDF output into PowerPoint.
    - This step has been integrated into the `02_run_sashimi_from_manifest.sh`, while this `03_build_sashimi_ppt_from_manifest.py` need to exist in the same folder that the program can read the script.
- `02-downstream-analysis/04_make_relative_frequency_logos_only.R`: creates relative nucleotide-frequency logos where supplied.



