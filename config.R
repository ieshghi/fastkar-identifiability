## Run configuration. Edit this file per cohort; every script sources it.

## -----------------------------------------------------------------------------
## Input and output
## -----------------------------------------------------------------------------

## Pairs table (.rds or tab-separated text): one row per sample, with columns
## `pair` (sample id) and `complex` (path to a gGnome graph .rds whose
## $meta$events has type, footprint and ev.id, i.e. gGnome::events() output).
pairs_file = '/path/to/pairs.rds'
cohort = 'COHORT'

## Everything the run writes goes under out_dir/run_tag: event table, job plan,
## bundles, Flow jobs (flow/), collected tables and plots (plots/).
out_dir = '/path/to/output'
run_tag = 'local_swap_v1'

## -----------------------------------------------------------------------------
## Events
## -----------------------------------------------------------------------------

## Simple one-junction event types are skipped.
exclude_types = c('dup', 'del', 'tra', 'inv')
## Set to a small integer for a test run on that many random events.
event_subset_n = Inf
## Per-event seeds are derived from this and the event id.
sample_seed = 1L
## Cores for reading the graphs when building the event table.
graph_cores = 8L

## -----------------------------------------------------------------------------
## Analysis settings, passed unchanged to fastKar:::run_graph_local_identifiability()
## (identical to the CCLE run in the manuscript, ccle_v3_local_swap)
## -----------------------------------------------------------------------------

settings = list(
    n_neighbors = 99,           # candidates = truth + up to 99 others
    n_walkers = 200,            # local.sampling() chains from the truth
    n_swaps = 3,                # swaps per chain; every visited karyotype is a candidate
    n_obs = 100,                # noisy observations of the truth per assay
    graph_mode = 'full',
    modalities = c('hic', 'longread'),
    longread_lengths = c(1e5, 1e7),
    depth = 10,
    theta = 2,
    hic_binsperft = 300,
    footprint_buffer_fraction = 0.2,
    footprint_buffer_min = 1e4,
    edit_mode = 'fragment',
    longread_background = 1e-10)

## -----------------------------------------------------------------------------
## Cluster jobs (Flow on Slurm)
## -----------------------------------------------------------------------------

## Cores per job (walkers and edit distances; scoring is single-core).
mc_cores = 4L
flow_update_cores = 4L
## Per-type runtime and memory from the CCLE run (4 cores per job):
## resources/ccle_event_resources.tsv. Types not in it use the 'ALL' row.
## Memory: smallest tier >= mem_scale x p90 RSS + mc_cores x fork_gb + mem_overhead_gb.
## Time: p90 minutes x time_scale. Events are packed into jobs of up to
## target_job_h estimated hours per memory tier; the Slurm limit is
## max(time_safety x estimate, estimate + time_pad_h). Failed jobs are rerun
## one memory tier up with twice the time by 02_submit.R.
mem_tiers_gb = c(16, 32, 64, 128, 256)
mem_scale = 1.5
fork_gb = 2
mem_overhead_gb = 4
time_scale = 2
target_job_h = 2
time_safety = 3
time_pad_h = 2

## -----------------------------------------------------------------------------
## Derived paths (no need to edit)
## -----------------------------------------------------------------------------

run_dir = file.path(path.expand(out_dir), run_tag)
flow_root = file.path(run_dir, 'flow')
plot_dir = file.path(run_dir, 'plots')
