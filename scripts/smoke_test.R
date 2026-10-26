## Run a few cheap events locally (no Slurm) to check the install and inputs
## before submitting.
##
##   source env.sh && Rscript scripts/smoke_test.R [--config=<file>] [--n=2] [--cores=4]
##
## Picks n events of the fastest types (CCLE estimates) from <run_dir>/events.rds
## (needs scripts/01_events.R), runs R/worker.R on them in <run_dir>/smoke/, and
## prints P(MLE = truth) per assay. A few minutes per event at 4 cores.
suppressPackageStartupMessages(library(data.table))
repo = normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))), '..'))
source(file.path(repo, 'R', 'common.R'))
load_config(repo)
arg_value = function(name, default) {
    a = grep(paste0('^--', name, '='), commandArgs(trailingOnly = TRUE), value = TRUE)
    if (length(a)) as.integer(sub('.*=', '', a[1])) else default
}
n = arg_value('n', 2L)
cores = arg_value('cores', mc_cores)

events = readRDS(run_file('events.rds'))
res = fread(file.path(repo, 'resources', 'ccle_event_resources.tsv'))
events[res, on = 'type', p90_min := i.p90_min]
set.seed(1)
events = events[order(p90_min, sample(.N))][seq_len(min(n, .N))]
message('Smoke test on: ', paste(events$event_id, ' (', events$type, ')', collapse = ', ', sep = ''))

smoke_dir = run_file('smoke')
unlink(smoke_dir, recursive = TRUE)
dir.create(smoke_dir)
saveRDS(events[, .(pair, event_id, cohort, type, ev.id, event.iid, identifiability_event_index, gg, footprint, seed)],
        file.path(smoke_dir, 'bundle.rds'))
saveRDS(settings, file.path(smoke_dir, 'settings.rds'))
old_wd = setwd(smoke_dir)
status = system2('Rscript', c('--vanilla', file.path(repo, 'R', 'worker.R'), 'smoke', 'bundle.rds', 'settings.rds', cores))
setwd(old_wd)
if (status != 0) stop('worker failed (exit ', status, '); see the output above and ', smoke_dir, '/events/*.failed')

picks = readRDS(file.path(smoke_dir, 'picks.rds'))
print(picks[, .(n_cand = n_cand[1], unique_walk = unique_walk[1], p_pick = mean(p_pick)),
            by = .(event_id, type, modality, readL)])
message('Smoke test passed. Outputs in ', smoke_dir)
