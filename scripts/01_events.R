## Build the event table from the pairs table's graphs.
##
##   source env.sh && Rscript scripts/01_events.R [--config=<file>] [--rebuild]
##
## One row per event: every event in each graph's $meta$events except
## exclude_types. event_id = <pair>_event_<index>, index running over the kept
## events in pairs-table order; each event gets a seed derived from sample_seed
## and its id, so its result does not depend on how events are bundled.
## With events_file set (config.R), only the listed (pair, ev.id) events are
## kept. Ids and seeds are still assigned over all events, so each listed event
## gets the same id, seed and result as in a full run; listed events of an
## excluded type are kept, with id <pair>_event_x<event.iid>. Extra columns of
## the events table go to <run_dir>/event_annotations.rds (and .tsv), which
## 03_analyze.R joins onto its tables.
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

all_events[, excluded := type %in% exclude_types]
all_events[!(excluded), identifiability_event_index := seq_len(.N)]
all_events[, event_id := ifelse(excluded, paste0(pair, '_event_x', event.iid),
                                paste0(pair, '_event_', identifiability_event_index))]
annotations = NULL
if (!exists('events_file')) events_file = NULL
if (is.null(events_file)) {
    events = all_events[!(excluded)]
    message('Dropped ', nrow(all_events) - nrow(events), ' events of type ', paste(exclude_types, collapse = '/'))
} else {
    wanted = read_event_list(events_file)
    key = function(p, e) paste(p, as.character(e), sep = '\r')
    all_key = key(all_events$pair, all_events$ev.id)
    dup = unique(all_events[duplicated(all_key) | duplicated(all_key, fromLast = TRUE)]$pair)
    if (length(dup)) stop('ev.id is not unique within the graph(s) of ', paste(head(dup, 5), collapse = ', '))
    hit = match(key(wanted$pair, wanted$ev.id), all_key)
    if (anyNA(hit)) {
        miss = wanted[is.na(hit)]
        miss[, why := ifelse(pair %in% pairs$pair, 'ev.id not in the graph\'s $meta$events', 'pair not in the pairs table')]
        print(head(miss[, .(pair, ev.id, why)], 20))
        stop(nrow(miss), ' of ', nrow(wanted), ' events in ', events_file, ' were not found (first 20 above)')
    }
    events = all_events[hit]
    message('Kept the ', nrow(events), ' events listed in ', events_file)
    if (any(events$excluded)) message(sum(events$excluded), ' of them are of an excluded type (',
                                      paste(unique(events[(excluded)]$type), collapse = '/'), '); running them anyway')
    extra = setdiff(names(wanted), c('pair', 'ev.id'))
    if (length(extra)) annotations = cbind(events[, .(event_id)], wanted[, ..extra])
}
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
ann_file = run_file('event_annotations.rds')
unlink(c(ann_file, run_file('event_annotations.tsv')))
if (!is.null(annotations)) {
    annotations = annotations[event_id %in% events$event_id]
    saveRDS(annotations, ann_file)
    fwrite(annotations, run_file('event_annotations.tsv'), sep = '\t')
    message('Wrote ', ann_file, ': columns ', paste(names(annotations)[-1], collapse = ', '))
}
message('Wrote ', out, ': ', nrow(events), ' events')
print(events[, .N, keyby = type])
