# No-op Flow profile. Flow (when run from Rscript) needs R_FLOW_PROFILE set,
# or it marks jobs "profile not ready" and writes empty cmd.sh files.
# scripts/02_submit.R uses this file unless R_FLOW_PROFILE is already set.
