## Plan, submit and rerun the Flow jobs (one job = a bundle of events run
## serially by R/worker.R).
##
##   source env.sh && Rscript scripts/02_submit.R [--config=<file>] [--dry-run]
##
## The first run packs events into jobs by estimated runtime within memory tiers
## (per-type estimates from resources/ccle_event_resources.tsv, see config.R),
## writes the bundles and settings, and submits. Rerunning reuses the saved plan,
## never resubmits jobs Slurm still holds, and resubmits failed jobs one memory
## tier up with twice the time (jobs that died without Flow noticing, e.g. Slurm
## timeouts, count as failed); finished events are kept (checkpoints), so only
## unfinished ones rerun. --dry-run plans and builds the Flow job, submits nothing.
##
## Needs scripts/01_events.R first. Collect and plot with scripts/03_analyze.R.
suppressPackageStartupMessages({
    library(data.table)
    library(Flow)
})
repo = normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))), '..'))
source(file.path(repo, 'R', 'common.R'))
load_config(repo)
setDTthreads(1)
submit_now = !has_flag('--dry-run')

module_dir = file.path(repo, 'flow', 'fastkar_local_identifiability')
## Flow needs a profile under Rscript, or it writes empty job commands.
if (!nzchar(Sys.getenv('R_FLOW_PROFILE'))) Sys.setenv(R_FLOW_PROFILE = file.path(repo, 'flow', 'flow_profile.sh'))
bundle_dir = run_file('bundles')
dir.create(bundle_dir, showWarnings = FALSE)
dir.create(flow_root, showWarnings = FALSE)

## Flow task: module path (this repo) plus inputs/outputs.
flow_task = run_file('fastkar_local_identifiability.task')
writeLines(c(paste0(module_dir, '/'),
             'input\tbatch_id\tbatch_id\tvalue',
             'input\tevent_bundle\tevent_bundle\tpath',
             'input\tsettings\tsettings\tpath',
             'input\tmc_cores\tmc_cores\tvalue',
             'output\tpicks_rds\tpicks.rds$',
             'output\tmeta_rds\tmeta.rds$',
             'output\tresult_rds\tresult.rds$',
             'output\tconfig_rds\tconfig.rds$'), flow_task)

## -----------------------------------------------------------------------------
## Plan: bundles, settings and Flow inputs (reused if already saved)
## -----------------------------------------------------------------------------

pick_tier = function(gb) vapply(gb, function(x) {
    t = mem_tiers_gb[mem_tiers_gb >= x]
    if (length(t)) min(t) else max(mem_tiers_gb)
}, numeric(1))

## First-fit decreasing on est_min.
pack_events = function(est_min, cap_min) {
    o = order(-est_min)
    load = numeric()
    bin = integer(length(est_min))
    for (i in o) {
        j = which(load + est_min[i] <= cap_min)[1]
        if (is.na(j)) { load = c(load, est_min[i]); j = length(load) }
        else load[j] = load[j] + est_min[i]
        bin[i] = j
    }
    bin
}

plan_files = run_file(c('plan_events.rds', 'plan_jobs.rds', 'plan_bundle_inputs.rds'))
settings_file = run_file('settings.rds')
if (all(file.exists(plan_files))) {
    message('Reusing the saved plan in ', run_dir)
    events = readRDS(plan_files[1])
    jobs = readRDS(plan_files[2])
    bundle_inputs = readRDS(plan_files[3])
} else {
    events_file = run_file('events.rds')
    if (!file.exists(events_file)) stop('no event table; run scripts/01_events.R first')
    events = readRDS(events_file)

    res = fread(file.path(repo, 'resources', 'ccle_event_resources.tsv'))
    events[res, on = 'type', `:=`(p90_min = i.p90_min, p90_rss_gb = i.p90_rss_gb)]
    unknown = unique(events[is.na(p90_min)]$type)
    if (length(unknown)) message('No CCLE estimate for type(s) ', paste(unknown, collapse = ', '), '; using all events')
    events[is.na(p90_min), `:=`(p90_min = res[type == 'ALL']$p90_min, p90_rss_gb = res[type == 'ALL']$p90_rss_gb)]
    events[, est_min := pmax(time_scale * p90_min, 0.5)]
    events[, mem_gb := pick_tier(mem_scale * p90_rss_gb + mc_cores * fork_gb + mem_overhead_gb)]
    events[, bin := pack_events(est_min, target_job_h * 60), by = mem_gb]
    events[, batch_id := sprintf('m%03d_%04d', as.integer(mem_gb), bin)]
    jobs = events[, .(n_events = .N, est_h = sum(est_min) / 60, mem_gb = mem_gb[1]), keyby = batch_id]
    jobs[, time_h := ceiling(pmax(time_safety * est_h, est_h + time_pad_h))]
    jobs[, time := sprintf('%d-%02d', time_h %/% 24, time_h %% 24)]

    saveRDS(settings, settings_file)
    bundle_inputs = rbindlist(lapply(split(events, by = 'batch_id'), function(dt) {
        bundle_path = file.path(bundle_dir, paste0(dt$batch_id[1], '.rds'))
        saveRDS(dt[, .(pair, event_id, cohort, type, ev.id, event.iid, identifiability_event_index,
                       gg, footprint, seed)], bundle_path)
        data.table(batch_id = dt$batch_id[1], event_bundle = bundle_path,
                   settings = settings_file, mc_cores = mc_cores)
    }))
    setkey(bundle_inputs, batch_id)
    saveRDS(events, plan_files[1])
    saveRDS(jobs, plan_files[2])
    saveRDS(bundle_inputs, plan_files[3])
}
message(nrow(events), ' events in ', nrow(jobs), ' jobs; est. ', round(sum(jobs$est_h)),
        ' job-hours at ', mc_cores, ' cores')
print(jobs[, .(jobs = .N, events = sum(n_events), est_job_h = round(sum(est_h)),
               max_est_h = round(max(est_h), 1), max_time_h = max(time_h)), keyby = mem_gb])

## -----------------------------------------------------------------------------
## Create/update/submit the Flow job
## -----------------------------------------------------------------------------

# Mem/time raised by earlier failed-job reruns. Read before Flow::Job(), which
# overwrites the cached job object.
jb_file = file.path(flow_root, 'fastkar_local_identifiability.rds')
prev = if (file.exists(jb_file)) readRDS(jb_file)@runinfo[, .(batch_id, mem, time)] else NULL

jb = Flow::Job(flow_task, bundle_inputs, rootdir = flow_root, mem = min(mem_tiers_gb),
               time = '0-04', update_cores = flow_update_cores, cores = mc_cores)
jb@runinfo[jobs$batch_id, mem := jobs$mem_gb]
jb@runinfo[jobs$batch_id, time := jobs$time]
if (!is.null(prev)) {
    prev = prev[batch_id %in% ids(jb)]
    jb@runinfo[prev$batch_id, `:=`(mem = prev$mem, time = prev$time)]
}
Flow::update(jb)

# Jobs Slurm still has (pending or running), by job name. Flow calls a queued
# job 'ready', so these are never submitted or rerun again.
live_ids = sub('^job\\.fastkar_local_identifiability\\.', '',
               system2('squeue', c('-u', Sys.getenv('USER'), '-h', '-o', '%j'), stdout = TRUE))
live_ids = intersect(live_ids, ids(jb))
if (length(live_ids)) message(length(live_ids), ' jobs are queued or running in Slurm')

# Longest jobs first, so the slowest events are not left for the end.
if (submit_now) {
    ready = setdiff(intersect(jobs[order(-est_h)]$batch_id, ids(jb['ready'])), live_ids)
    if (length(ready)) {
        message('Submitting ', length(ready), ' ready jobs')
        Flow::srun(jb[ready, id = TRUE])
    } else {
        message('No ready jobs to submit')
    }
    Flow::update(jb)
}
saveRDS(jb, jb_file)

## Failed jobs: rerun at the next memory tier (doubled past the top tier) and
## twice the time. Bundle outputs and Slurm files are cleared; the events/
## checkpoints are kept, so finished events are not redone. A job Flow still
## calls running whose Slurm job is gone (timeout, node failure, cancelled)
## counts as failed.
failed_ids = setdiff(ids(jb)[grep('failed', status(jb))], live_ids)
stale = setdiff(ids(jb)[grep('running', status(jb))], live_ids)
if (length(stale)) message(length(stale), ' running jobs have no live Slurm job; treating as failed: ',
                           paste(stale, collapse = ', '))
failed_ids = c(failed_ids, stale)
if (length(failed_ids) && !submit_now) {
    message('--dry-run; would rerun: ', paste(failed_ids, collapse = ', '))
} else if (length(failed_ids)) {
    ri = jb@runinfo[failed_ids]
    new_mem = vapply(ri$mem, function(m) {
        t = mem_tiers_gb[mem_tiers_gb > m]
        if (length(t)) min(t) else 2 * m
    }, numeric(1))
    old_h = vapply(strsplit(as.character(ri$time), '-', fixed = TRUE),
                   function(x) 24 * as.numeric(x[1]) + as.numeric(x[2]), numeric(1))
    new_h = 2 * old_h
    message('Rerunning ', length(failed_ids), ' failed jobs at higher mem/time')
    jb@runinfo[failed_ids, mem := new_mem]
    jb@runinfo[failed_ids, time := sprintf('%d-%02d', new_h %/% 24, new_h %% 24)]
    for (d in outdir(jb[failed_ids, id = TRUE])) {
        unlink(file.path(d, c('picks.rds', 'meta.rds', 'result.rds', 'config.rds',
                              'started', 'slurm.jobid')), force = TRUE)
        unlink(Sys.glob(file.path(d, c('*.cmd.sh', '*.bsub.out', '*.bsub.err', '*.bsub.report'))), force = TRUE)
    }
    Flow::update(jb, check.inputs = FALSE, mc.cores = flow_update_cores, cache.object = TRUE)
    Flow::srun(jb[failed_ids, id = TRUE])
    Flow::update(jb, check.inputs = FALSE, mc.cores = flow_update_cores, cache.object = TRUE)
    saveRDS(jb, jb_file)
}

message('Status:')
print(table(status(jb)))
message('Collect and plot with scripts/03_analyze.R')
