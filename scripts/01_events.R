## Build the event table from the pairs table's graphs.
##
##   source env.sh && Rscript scripts/01_events.R [--config=<file>] [--rebuild]
##
## One row per event: every event in each graph's $meta$events except
## exclude_types. event_id = <pair>_event_<index>, index running over the kept
## events in pairs-table order; each event gets a seed derived from sample_seed
## and its id, so its result does not depend on how events are bundled.
## Writes <run_dir>/events.rds (and events.tsv). An existing table is kept
## unless --rebuild is given: event ids and seeds must not change under a live run.
suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
})
repo = normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))), '..'))
source(file.path(repo, 'R', 'common.R'))
load_config(repo)
setDTthreads(1)

out = run_file('events.rds')
if (file.exists(out) && !has_flag('--rebuild')) {
    message(out, ' exists; keeping it (use --rebuild to remake it)')
    quit(save = 'no')
}

pairs = read_pairs(pairs_file)
missing = pairs[!file.exists(path.expand(complex))]
if (nrow(missing)) stop(nrow(missing), ' graph files do not exist, e.g. ', missing$complex[1])
message('Reading events from ', nrow(pairs), ' graphs')
evs = mclapply(seq_len(nrow(pairs)), function(i) {
    gg = readRDS(path.expand(pairs$complex[i]))
    ev = gg$meta$events
    if (is.null(ev) || !nrow(ev)) return(NULL)
    miss = setdiff(c('type', 'footprint', 'ev.id'), names(ev))
    if (length(miss)) stop(pairs$pair[i], ': $meta$events lacks ', paste(miss, collapse = ', '))
    ev = as.data.table(ev)[, .(type, footprint, ev.id)]
    ev[, pair := pairs$pair[i]]
    ev
}, mc.cores = graph_cores)
bad = vapply(evs, inherits, logical(1), 'try-error')
if (any(bad)) stop('reading graphs failed: ', as.character(evs[bad][[1]]))
all_events = rbindlist(evs, fill = TRUE)[, event.iid := .I][, cohort := cohort]
message(nrow(all_events), ' events in ', uniqueN(all_events$pair), ' samples (',
        nrow(pairs) - uniqueN(all_events$pair), ' samples have none)')

events = all_events[!(type %in% exclude_types)]
message('Dropped ', nrow(all_events) - nrow(events), ' events of type ', paste(exclude_types, collapse = '/'))
events[, identifiability_event_index := .I]
events[, event_id := paste0(pair, '_event_', identifiability_event_index)]
events[, gg := pairs[match(events$pair, pairs$pair), complex]]
if (is.finite(event_subset_n)) {
    set.seed(sample_seed)
    events = events[sort(sample(.N, min(.N, event_subset_n)))]
    message('Test run: kept ', nrow(events), ' random events')
}
events[, seed := vapply(event_id, function(e) strtoi(substr(digest::digest(
    paste(sample_seed, e, 'local'), algo = 'xxhash32', serialize = FALSE), 1, 7), 16L), integer(1))]
events = events[, .(pair, event_id, cohort, type, ev.id, event.iid, identifiability_event_index, gg, footprint, seed)]

saveRDS(events, out)
fwrite(events, run_file('events.tsv'), sep = '\t')
message('Wrote ', out, ': ', nrow(events), ' events')
print(events[, .N, keyby = type])
