# Software environment

The package lists in this directory were extracted from the supplied scripts. They identify dependencies but do not yet reproduce the original software environment because exact versions were not supplied.

Before publication:

1. Restore or load the final analysis environment.
2. Run `sessionInfo()` after loading the packages used by each analysis.
3. Create and commit an `renv.lock` file if the final scripts can be restored with `renv`.
4. Pin Python package versions in `python-requirements.txt`.
5. Record the container digest or immutable version for rMATS and other command-line tools.

Do not generate a lockfile from an unrelated modern environment; it should represent the environment used to produce or validate the reported results.
