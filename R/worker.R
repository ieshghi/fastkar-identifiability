## Local identifiability for a bundle of events (fastKar:::run_graph_local_identifiability()).
##
##   Rscript R/worker.R <batch_id> <event_bundle.rds> <settings.rds> <mc_cores>
##
## Run by the Flow module (flow/fastkar_local_identifiability/run.sh) in its job
## directory, or by hand in any directory (scripts/smoke_test.R does this).
## Each finished event is checkpointed to events/<event_id>.rds, so a rerun
## (e.g. after an OOM kill) skips it. An event that errors writes
## events/<event_id>.failed; the job then fails after trying the rest of the
## bundle, and no bundle outputs are written.
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop('usage: worker.R <batch_id> <event_bundle> <settings> <mc_cores>')
batch_id = args[1]
event_bundle = path.expand(args[2])
settings_file = path.expand(args[3])
mc.cores = as.integer(args[4])
slurm_cpus = suppressWarnings(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')))
if (!is.na(slurm_cpus) && mc.cores > slurm_cpus) {
    message('mc_cores = ', mc.cores, ' but SLURM_CPUS_PER_TASK = ', slurm_cpus, '; using ', slurm_cpus)
    mc.cores = slurm_cpus
}
if (!is.finite(mc.cores) || mc.cores < 1) stop('mc_cores must be positive')

suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
    library(fastKar)
    library(gUtils)
    # gGnome after fastKar, as in interactive use.
    library(gGnome)
})
setDTthreads(1)

## Some gUtils versions call aggregate(formula = ...) in grl.in(), which newer R
## rejects; patch it in place if so (no-op otherwise).
patch_gutils_grl_in = function() {
    ns = asNamespace('gUtils')
    if (!exists('grl.in', envir = ns, inherits = FALSE)) return(invisible(FALSE))
    old = get('grl.in', envir = ns)
    b = paste(deparse(body(old)), collapse = '\n')
    b2 = sub('stats::aggregate\\(formula = subject.id ~ grl.id,\\s*data = m, FUN = function\\(x\\) numwin - length\\(setdiff\\(1:numwin,\\s*x\\)\\)\\)',
             'stats::aggregate(subject.id ~ grl.id, data = m, FUN = function(x) numwin - length(setdiff(1:numwin, x)))', b)
    if (identical(b, b2)) return(invisible(FALSE))
    new = old
    body(new) = parse(text = b2)[[1]]
    environment(new) = ns
    unlockBinding('grl.in', ns); assign('grl.in', new, envir = ns); lockBinding('grl.in', ns)
    if ('package:gUtils' %in% search()) {
        pkg_env = as.environment('package:gUtils')
        unlockBinding('grl.in', pkg_env); assign('grl.in', new, envir = pkg_env); lockBinding('grl.in', pkg_env)
    }
    message('Patched gUtils::grl.in aggregate(formula=...) compatibility issue')
    invisible(TRUE)
}
patch_gutils_grl_in()

## Primary chromosomes (1-22, X, Y, with or without 'chr') of the graph's own
## genome: footprints are clipped to these, so any build and naming works.
graph_genome = function(gg) {
    gr = gg$nodes$gr
    sl = GenomeInfoDb::seqlengths(gr)
    std = names(sl)[grepl('^(chr)?([0-9]+|X|Y)$', names(sl))]
    if (!length(std)) stop('graph has no primary chromosomes (1-22, X, Y)')
    if (all(!is.na(sl[std]))) return(GenomicRanges::GRanges(std, IRanges::IRanges(1, sl[std])))
    streduce(gr[as.character(GenomicRanges::seqnames(gr)) %in% std])
}

settings = readRDS(settings_file)
bundle = as.data.table(readRDS(event_bundle))
if (!nrow(bundle)) stop('event_bundle has zero rows')
required_cols = c('pair', 'event_id', 'cohort', 'type', 'gg', 'footprint', 'seed')
missing_cols = setdiff(required_cols, names(bundle))
if (length(missing_cols)) stop('event_bundle missing columns: ', paste(missing_cols, collapse = ', '))
if (anyDuplicated(bundle$event_id)) stop('event_bundle has duplicated event_id rows')

event_dir = 'events'
dir.create(event_dir, showWarnings = FALSE)
event_file = function(e, ext = '.rds') file.path(event_dir, paste0(e, ext))
done = file.exists(event_file(bundle$event_id))
message('Running local identifiability batch ', batch_id, ': ', nrow(bundle), ' events, ',
        sum(done), ' already done (checkpoints)')

# L'Ecuyer streams make the forked local.sampling() walkers follow from each
# event's seed (for a fixed mc_cores).
RNGkind("L'Ecuyer-CMRG")
# Events run one at a time and only their results are saved, so the parent
# heap (copied into every forked worker) does not grow across the bundle.
for (k in which(!done)) {
    ev = bundle[k]
    unlink(event_file(ev$event_id, '.failed'))
    r = tryCatch({
        gg = readRDS(normalizePath(path.expand(ev$gg), mustWork = TRUE))
        if (!is.null(gg$copy)) gg = gg$copy
        do.call(fastKar:::run_graph_local_identifiability, c(list(
            graph = gg,
            footprint = streduce(parse.gr(strsplit(ev$footprint, ';', fixed = TRUE)[[1]])),
            event = ev[, .(event_id, type, pair, cohort)],
            seed = ev$seed,
            wholegenome = graph_genome(gg),
            mc.cores = mc.cores), settings))
    }, error = function(e) e)
    if (inherits(r, 'error')) {
        message('[', ev$event_id, '] FAILED: ', conditionMessage(r))
        writeLines(conditionMessage(r), event_file(ev$event_id, '.failed'))
    } else {
        saveRDS(r, event_file(ev$event_id))
    }
    rm(r); invisible(gc())
}

failed = bundle$event_id[!file.exists(event_file(bundle$event_id))]
if (length(failed)) stop(length(failed), ' of ', nrow(bundle), ' events failed: ', paste(failed, collapse = ', '))

res = lapply(bundle$event_id, function(e) readRDS(event_file(e)))
names(res) = bundle$event_id
picks = rbindlist(lapply(res, `[[`, 'picks'), fill = TRUE)
meta = rbindlist(lapply(res, `[[`, 'meta'), fill = TRUE)
config = list(batch_id = batch_id, event_bundle = event_bundle, settings_file = settings_file,
              settings = settings, inputs = bundle, mc_cores = mc.cores,
              packages = sapply(c('fastKar', 'gGnome', 'gUtils', 'GxG'), function(p) as.character(packageVersion(p))))
saveRDS(res, 'result.rds')
saveRDS(meta, 'meta.rds')
saveRDS(config, 'config.rds')
# picks.rds last: its presence marks the job complete.
saveRDS(picks, 'picks.rds')
message('Local identifiability batch ', batch_id, ' complete')
