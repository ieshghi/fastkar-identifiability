## Install fastKar and its dependencies into the first library path
## (FASTKAR_R_LIB from env.sh). Run once per cluster:
##
##   source env.sh && Rscript install.R            # pinned versions below
##   source env.sh && Rscript install.R --flow     # also install Flow
##
## CRAN and Bioconductor dependencies are installed as needed; packages already
## present are not upgraded. The GitHub packages are pinned to the commits this
## pipeline was tested with (CCLE checkpoints reproduced exactly).

lib = .libPaths()[1]
if (!nzchar(Sys.getenv('FASTKAR_R_LIB')) || normalizePath(lib) != normalizePath(Sys.getenv('FASTKAR_R_LIB')))
    stop('source env.sh first, so packages go into FASTKAR_R_LIB and not your default library')
message('Installing into ', lib)

repos = c(CRAN = 'https://cloud.r-project.org')
for (p in c('remotes', 'BiocManager'))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib, repos = repos)
repos = BiocManager::repositories()

## Order matters: each needs the ones before it. gUtils, gGnome and GxG must be
## installed together (GxG needs a gUtils newer than some cluster installs).
pinned = c(
    'mskilab-org/gUtils@095377a',
    'mskilab-org/gTrack@2ec3da3',
    'mskilab-org/skitools@ba322dc',
    'mskilab-org/gGnome@803c279',
    'mskilab-org/GxG@a8d5827',
    'ieshghi/fastKar@1c5197b')
if ('--flow' %in% commandArgs(trailingOnly = TRUE)) pinned = c(pinned, 'mskilab-org/Flow@63ce03d')

for (r in pinned) {
    message('\n== ', r)
    remotes::install_github(r, lib = lib, repos = repos, dependencies = TRUE,
                            upgrade = 'never', force = TRUE)
}
## Used only by scripts/03_analyze.R for the figure layout.
for (p in c('ggforce', 'egg'))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib, repos = repos)

suppressPackageStartupMessages(library(fastKar))
message('\nfastKar ', packageVersion('fastKar'), ' loads from ', find.package('fastKar'))
