#!/bin/bash
# Flow module: local identifiability for one bundle of events.
# Args: batch_id event_bundle settings mc_cores (see R/worker.R).
MODULE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$MODULE_DIR/../.." && pwd)
source "$REPO/env.sh"
exec Rscript --vanilla "$REPO/R/worker.R" "$@"
