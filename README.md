# fastkar-identifiability

Pan-cancer **local karyotype identifiability** with [fastKar](https://github.com/ieshghi/fastKar),
as in Fig. 4 of *The identifiability frontier of cancer karyotype assembly*.

For each rearrangement event in a cohort's genome graphs: draw a random truth karyotype, collect up to
99 nearby karyotypes (1-3 local swaps at footprint nodes), simulate 100 noisy observations of the
truth per assay (Hi-C; long reads at 100 kb and 10 Mb; all at 10x), and record
**P(MLE = truth)**, the probability that each assay's maximum-likelihood karyotype is the truth
(ties broken uniformly; chance = 1 / number of candidates). Everything is simulated from the graphs:
no sequencing data are needed.

## Input

A pairs table (`.rds`, or tab-separated text) with one row per sample and columns

| column | content |
|---|---|
| `pair` | sample id |
| `complex` | path to a gGnome graph `.rds` with events annotated: `$meta$events` with `type`, `footprint`, `ev.id` (`gGnome::events()` output, as in JaBbA `complex.rds`) |

Any genome build and chromosome naming works: footprints are clipped to each graph's own
primary chromosomes (1-22, X, Y, with or without `chr`).

## Setup (once per cluster)

Needs R >= 4.3 with a C++ compiler, Slurm and [Flow](https://github.com/mskilab-org/Flow).

1. Edit the block at the top of `env.sh`: how to get R (e.g. `module load R/4.3.3`) and
   `FASTKAR_R_LIB`, the R library this pipeline installs into (default `~/R/fastkar-identifiability`,
   kept separate from your own library).
2. Install fastKar and its GitHub dependencies (gUtils, gTrack, skitools, gGnome, GxG) at the pinned
   commits in `install.R`; CRAN/Bioconductor dependencies come along:
   ```bash
   source env.sh && Rscript install.R          # add --flow to install Flow too
   ```

## Run

Edit `config.R` (pairs table, cohort name, output directory; analysis settings default to the
manuscript's CCLE run), then:

```bash
source env.sh
Rscript scripts/01_events.R        # event table from the graphs -> <out_dir>/<run_tag>/events.rds
Rscript scripts/smoke_test.R       # 2 cheap events locally, a few minutes: checks install and inputs
Rscript scripts/02_submit.R --dry-run   # plan jobs, build the Flow job, submit nothing
Rscript scripts/02_submit.R        # submit
Rscript scripts/02_submit.R        # later: rerun failed jobs one memory tier up (finished events are kept)
Rscript scripts/03_analyze.R       # any time: collect finished events, tables and plots
```

All scripts take `--config=<file>` to use another config file.

### Running a chosen set of events

To run only some events, e.g. those overlapping a gene set or called by another method, list them
in an events table and point `events_file` in `config.R` at it:

```r
events_file = '/path/to/my_events.tsv'   # NULL (default) = all events
```

The table (`.rds` or tab-separated text) has one row per event:

| column | content |
|---|---|
| `pair` | sample id, as in the pairs table |
| `ev.id` | the event's id in that sample's graph (`gg$meta$events$ev.id`) |
| any others | optional annotations, e.g. `group`, carried into the results |

```
pair	ev.id	group
G32831.LU99	5	focal
G32831.HCC38	17	focal
G32831.PA-TU-8988T	31	control
```

Every listed event must be in its graph's `$meta$events`; `01_events.R` stops and lists any that
are not, saying whether the `pair` is missing from the pairs table or the `ev.id` from the graph.
The rest of the pipeline is unchanged (`smoke_test.R`, `02_submit.R`, `03_analyze.R`).

- **Same ids, seeds and results as a full run.** Event ids and seeds are assigned over all events in
  the pairs table before the list is applied, so a listed event gets exactly the result it would get
  in a full run with the same pairs table and `sample_seed` (and can be compared or merged with one).
- **`exclude_types` does not apply to listed events**: a listed `dup`, `del`, `tra` or `inv` is run
  (with id `<pair>_event_x<n>`, as it has no id in a full run).
- **`event_subset_n`** still works and samples from the listed events.
- **Annotations.** Extra columns are saved as `event_annotations.rds/.tsv` and joined onto
  `event_summary.rds/.tsv` by `03_analyze.R`, e.g. to compare groups.
- **Use a new `run_tag` for each event set.** `01_events.R` keeps an existing event table and
  `02_submit.R` reuses a saved job plan, so changing the list under an existing `run_tag` does nothing
  unless you start a new run directory.

### Outputs (`<out_dir>/<run_tag>/`)

| file | content |
|---|---|
| `events.rds/.tsv` | events analysed (one row per event; `event_id = <pair>_event_<n>`, per-event `seed`) |
| `event_annotations.rds/.tsv` | extra columns of the events table, by `event_id` (only with `events_file`) |
| `plan_*.rds`, `bundles/`, `settings.rds` | job plan, per-job event bundles, analysis settings |
| `flow/fastkar_local_identifiability/<job>/events/<event_id>.rds` | per-event result: `picks`, `meta`, `edit`, `swaps`, `nll` (NLL matrices, truth = column 1) |
| `picks.rds` | one row per event x assay x observation (truth rank, ties, `p_pick`) |
| `event_summary.rds/.tsv` | one row per event x assay: `p_pick` = P(MLE = truth), `chance`, `n_cand`, `unique_walk`, plus any event annotations |
| `plots/identifiability_by_type.pdf` | P(MLE = truth) by event type and assay |
| `plots/identifiability_figure.pdf` | manuscript Fig. 4 layout |

Events whose footprint admits only one karyotype (`unique_walk`) count as P = 1 for every assay.

## Cluster resources

Jobs run bundles of events serially at `mc_cores` (default 4) cores. Memory tier and runtime per event
are estimated from its type, using the CCLE run (`resources/ccle_event_resources.tsv`: per-type
median / 90th percentile / max minutes and parent-process RSS at 4 cores, 2233 events); see the
cluster block in `config.R`. Underestimates are handled by rerunning `02_submit.R`, which resubmits
failed jobs one memory tier up with twice the time; finished events are checkpointed and not redone.
On CCLE (2239 events) the plan is 242 jobs, ~500 estimated job-hours (the estimate is conservative).

## Layout

```
config.R                 run configuration (edit per cohort; optional events_file)
env.sh                   cluster environment (edit per cluster)
install.R                pinned package install
R/worker.R               runs one bundle of events (called by the Flow module or smoke_test.R)
R/common.R               helpers shared by the scripts
scripts/01_events.R      event table
scripts/02_submit.R      plan, submit, rerun
scripts/03_analyze.R     collect, summarise, plot
scripts/smoke_test.R     local test on a few events
flow/                    Flow module (run.sh -> R/worker.R) and no-op Flow profile
resources/               CCLE per-type runtime/memory
```

## Reproducibility

The pinned fastKar (`1c5197b`) differs from the one used for the manuscript's CCLE run
(`ccle_v3_local_swap`, fastKar `91f1ff3`) in two ways that change Hi-C NLLs:
- each pixel's Hi-C NLL is weighted by its area (width_i * width_j / pix.size^2; the tiles are cut at node
  boundaries, so pixels differ in size; since `bd9941c`), where `91f1ff3` scored every pixel equally;
- the forward model keeps the full length of nodes that the event's target region only partly covers
  (since `1c5197b`). Before, the uncovered part of such a node was dropped from the walk, so tiles on either
  side of a gap in the target sat next to each other. Nearly every event has such a node at the edge of its
  buffered footprint; about a fifth have a gap inside one node, where the effect is large. On a 42-event
  sample, 12 events changed P(MLE = truth) by more than 0.05 (mean 0.79 -> 0.78). With the default settings
this pipeline gives the same 2239 events, ids and seeds as that run, and per event the same candidates,
edit distances, noisy observations and long-read results (checked on a production checkpoint); only the
Hi-C NLLs, and so the Hi-C picks, can differ. `03_analyze.R` on the production outputs gives the
published numbers (mean P(MLE = truth): Hi-C 0.86, 10 Mb 0.73, 100 kb 0.35). Install `ieshghi/fastKar@91f1ff3`
instead to reproduce that run exactly. Results depend on
`mc_cores` (each event's forked sampling walkers draw from L'Ecuyer streams split by core).
