#!/bin/bash
#SBATCH -p tsl-medium
#SBATCH -t 1-00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -J guy11_bootstrap
#SBATCH -o results/hpc_logs/guy11_bootstrap.%j.out
#SBATCH -e results/hpc_logs/guy11_bootstrap.%j.err

# Guy11 within-experiment time-course bootstrap-t (10,000 iterations),
# then between-experiment kinetics correlation. See
# analysis/scripts/phosphoproteome_diffexp_bootstrap.r for the analysis
# rationale, and run_guy11_bootstrap_diffexp.r for what this actually runs.
#
# BEFORE SUBMITTING -- things I can't verify from here, please check:
#
#   1. Partition/time limit: -p tsl-medium and -t 1-00:00 (1 day) are a
#      guess with generous buffer -- estimated actual runtime is roughly
#      2-4 hours with -c 8 (16,155 shared Guy11 phosphosites x 5 timepoints
#      x 2 experiments x 10,000 bootstrap iterations each, and every one of
#      them runs the FULL bootstrap since the imputed input data has no
#      missing values to short-circuit on). Adjust -p/-t/-c to whatever
#      your site's actual queue limits and node core counts are.
#
#   2. R module: uncomment and fix the `source` line below to whatever
#      your site calls its R 4.3.x module (must be R 4.3 to match this
#      project's renv/library/R-4.3 -- a different R version's binary
#      package ABI won't load). Check with your site's software list /
#      sysadmin if you don't already know the package name.
#
#   3. renv::restore() below will very likely try to rebuild `pepdiff` from
#      https://github.com/TeamMacLean/pepdiff (recorded in renv.lock), and
#      that package's DESCRIPTION lists ARTool as a hard dependency it
#      doesn't actually need for what this script uses (only for pepdiff's
#      unused "ART" analysis method). On this laptop, ARTool's own
#      dependency chain (car -> pbkrtest -> doBy -> Deriv) failed to
#      compile against the local R/Rcpp headers. That may or may not
#      happen again on the HPC's toolchain -- if renv::restore() fails on
#      ARTool/Deriv/pepdiff, the fix that worked here was: download
#      https://github.com/TeamMacLean/pepdiff, delete the "ARTool," line
#      from its DESCRIPTION's Imports, then `R CMD INSTALL` that local copy
#      directly (bypassing renv) instead of relying on restore() for just
#      that one package.
#
# source R-4.3.0   # <-- fix this to your site's actual R module name

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p results/hpc_logs

echo "Host: $(hostname)  |  Job: $SLURM_JOB_ID  |  CPUs: $SLURM_CPUS_PER_TASK  |  Started: $(date)"

# Rebuild the renv library for this machine (see caveat #3 above).
Rscript -e 'renv::restore(prompt = FALSE)'

Rscript analysis/scripts/run_guy11_bootstrap_diffexp.r

echo "Finished: $(date)"
