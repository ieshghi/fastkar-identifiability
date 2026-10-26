# Cluster environment: an R >= 4.3 with a C++ compiler, and the R library that
# holds fastKar and its dependencies (filled by install.R). Edit the block
# below per cluster. Source this before running any script:
#
#   source env.sh && Rscript scripts/01_events.R
#
# The Flow jobs source it too (flow/fastkar_local_identifiability/run.sh).

## ---- edit per cluster ----------------------------------------------------------
# module load R/4.3.3
# R library for this pipeline (install.R installs into the first entry).
export FASTKAR_R_LIB=${FASTKAR_R_LIB:-$HOME/R/fastkar-identifiability}
## -----------------------------------------------------------------------------------

mkdir -p "$FASTKAR_R_LIB"
export R_LIBS="$FASTKAR_R_LIB${R_LIBS:+:$R_LIBS}"
# Keep the user's ~/.Rprofile out of every R process run by the pipeline.
export R_PROFILE_USER=/dev/null
# One thread per process: parallelism comes from forked workers (mc_cores).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
