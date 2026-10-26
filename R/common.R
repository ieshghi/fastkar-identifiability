## Shared by the scripts in scripts/. Sourced after `repo` is set.

## Source the run configuration into the global environment: config.R in the
## repo, or the file given as --config=<path>.
load_config = function(repo) {
    arg = grep('^--config=', commandArgs(trailingOnly = TRUE), value = TRUE)
    f = if (length(arg)) sub('^--config=', '', arg[1]) else file.path(repo, 'config.R')
    message('Config: ', normalizePath(f))
    sys.source(f, envir = globalenv())
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    invisible(f)
}

## Pairs table from .rds or delimited text; needs `pair` and `complex`.
read_pairs = function(f) {
    pairs = if (grepl('\\.rds$', f, ignore.case = TRUE)) readRDS(f) else data.table::fread(f)
    pairs = data.table::as.data.table(pairs)
    miss = setdiff(c('pair', 'complex'), names(pairs))
    if (length(miss)) stop('pairs table lacks column(s): ', paste(miss, collapse = ', '))
    if (anyDuplicated(pairs$pair)) stop('pairs table has duplicated `pair` values')
    pairs
}

## Run files, all under run_dir.
run_file = function(name) file.path(run_dir, name)

has_flag = function(flag) flag %in% commandArgs(trailingOnly = TRUE)
