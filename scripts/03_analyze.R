## Collect and plot the run. Safe to run at any time: it reads every finished
## event checkpoint in the Flow job directories (<job dir>/events/<event_id>.rds),
## so events from running or failed jobs are included too.
##
##   source env.sh && Rscript scripts/03_analyze.R [--config=<file>]
##
## Outputs (in run_dir):
##   picks.rds          one row per event x assay x observation: truth rank, ties, MLE pick
##   event_summary.rds  one row per event x assay (means over observations); also .tsv
## Plots (in plot_dir):
##   identifiability_by_type.pdf   P(MLE = truth) per event type and assay
##   identifiability_figure.pdf    manuscript Fig. 4 layout (A simple types, B complex
##                                 types, C per-sample mean, D fraction identifiable)
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggforce)
})
repo = normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))), '..'))
source(file.path(repo, 'R', 'common.R'))
load_config(repo)
setDTthreads(4)
# egg::ggarrange opens a device; keep it from writing Rplots.pdf in the working directory.
pdf(NULL)

event_dirs = Sys.glob(file.path(flow_root, 'fastkar_local_identifiability', '*', 'events'))
files = list.files(event_dirs, pattern = '\\.rds$', full.names = TRUE)
files = files[!duplicated(basename(files))]
failed = setdiff(sub('\\.failed$', '', list.files(event_dirs, pattern = '\\.failed$')),
                 sub('\\.rds$', '', basename(files)))
n_planned = if (file.exists(f <- run_file('plan_events.rds'))) nrow(readRDS(f)) else NA
message(length(files), ' of ', n_planned, ' events finished, ', length(failed), ' failed (not yet rerun)')
if (!length(files)) stop('nothing to collect yet')
res = parallel::mclapply(files, function(f) {
    r = readRDS(f)
    list(picks = r$picks, meta = r$meta)
}, mc.cores = 8)
picks = rbindlist(lapply(res, `[[`, 'picks'), fill = TRUE)
meta = rbindlist(lapply(res, `[[`, 'meta'), fill = TRUE)
picks[, assay := fcase(modality == 'hic', 'Hi-C',
                       readL == 1e5, 'LR 100 kb',
                       readL == 1e7, 'LR 10 Mb')]
picks[is.na(assay), assay := paste0('LR ', readL / 1e3, ' kb')]
picks[, assay := factor(assay, unique(c('LR 100 kb', 'LR 10 Mb', sort(unique(assay)), 'Hi-C')))]
picks[, assay := droplevels(assay)]

# p_pick: probability the MLE is the truth (ties broken uniformly).
# chance = 1 / n_cand is p_pick for an assay that sees nothing.
# unique_walk: the event admits only one karyotype inside the footprint, so
# p_pick = 1 for every assay (one row per assay instead of n_obs).
summ = picks[, .(p_pick = mean(p_pick), median_rank = median(as.numeric(truth_rank)),
                 chance = 1 / n_cand[1], n_obs = .N),
             by = .(event_id, pair, cohort, type, assay, n_cand, unique_walk, event_size, target_width)]
saveRDS(picks, run_file('picks.rds'))
saveRDS(summ, run_file('event_summary.rds'))
fwrite(summ, run_file('event_summary.tsv'), sep = '\t')

message('\nCandidate sets: ', sum(meta$unique_walk), ' events with a unique walk (p_pick = 1); median ',
        median(meta$n_cand), ' candidates')
message('\nPer-assay means over all events:')
print(summ[, .(events = .N, p_pick = mean(p_pick), frac_p_pick_gt_0.5 = mean(p_pick > 0.5),
               chance = mean(chance), median_rank = median(median_rank)), keyby = assay])

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
save_plot = function(p, name, width = 5, height = 3) {
    ggsave(file.path(plot_dir, paste0(name, '.pdf')), p, width = width, height = height)
    message('wrote ', file.path(plot_dir, name), '.pdf')
}
save_plot(ggplot(summ, aes(reorder(type, p_pick), p_pick, color = assay)) +
          geom_boxplot(aes(fill = assay), color = 'grey', alpha = 0.1) + geom_sina(size = 0.5) +
          labs(y = 'P(MLE = truth)', x = NULL, color = 'Assay') + theme_bw(),
          'identifiability_by_type', width = 10, height = 3)

## Manuscript Fig. 4 layout (needs the default assays: 100 kb, 10 Mb, Hi-C).
## qrpmin/qrppos/qrpmix are reciprocal duplications (r-dup); cpxdm is folded
## into dm. Types outside simple_types/complex_types appear only in the plot above.
s = copy(summ)
if (!all(c('LR 100 kb', 'LR 10 Mb', 'Hi-C') %in% s$assay)) {
    message('Skipping the figure layout: it needs the 100 kb, 10 Mb and Hi-C assays')
    quit(save = 'no')
}
s = s[assay %in% c('LR 100 kb', 'LR 10 Mb', 'Hi-C')]
s[type %in% c('qrpmin', 'qrppos', 'qrpmix'), type := 'r-dup']
s[type == 'cpxdm', type := 'dm']
simple_types = intersect(c('invdup', 'r-dup', 'tic', 'rigma', 'pyrgo'), s$type)
complex_types = intersect(c('chromoplexy', 'bfb', 'chromothripsis', 'dm', 'tyfonas'), s$type)
s[, assay_short := factor(sub('LR ', '', sub(' kb', 'kb', sub(' Mb', 'Mb', assay))),
                          c('100kb', '10Mb', 'Hi-C'))]
assay_cols = c('100kb' = '#1baf7a', '10Mb' = '#eb6834', 'Hi-C' = '#2a78d6')
ident_cutoff = 0.9
theme_fig4 = theme_minimal(base_size = 11) +
    theme(legend.position = 'none',
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = 'grey90', linewidth = 0.3),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.9),
          axis.ticks = element_line(color = 'grey30', linewidth = 0.3),
          axis.text = element_text(color = 'grey20'),
          strip.background = element_blank(),
          strip.text = element_text(size = 13, face = 'bold', color = 'grey10'),
          panel.spacing = unit(0, 'lines'))
type_panel = function(types, chance_label = FALSE) {
    d = s[type %in% types]
    n = d[assay == 'Hi-C', .N, by = type]
    d[n, label := paste0(i.type, ' (', i.N, ')'), on = 'type']
    d[, label := factor(label, n[match(types, type), paste0(type, ' (', N, ')')])]
    ch = d[, .(chance = mean(chance)), by = label]
    p = ggplot(d, aes(assay_short, p_pick, color = assay_short)) +
        geom_hline(data = ch, aes(yintercept = chance), linetype = 'dashed', color = 'grey35', linewidth = 0.4) +
        geom_sina(size = 0.5, alpha = 0.45, stroke = 0) +
        geom_boxplot(aes(fill = assay_short), outlier.shape = NA, alpha = 0.12, width = 0.55, linewidth = 0.45) +
        facet_wrap(~label, nrow = 1) +
        scale_color_manual(values = assay_cols) + scale_fill_manual(values = assay_cols) +
        scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = expansion(add = 0.02)) +
        labs(x = NULL, y = 'Identifiability\nP(MLE = truth)') + theme_fig4
    if (chance_label) {
        ## label the dashed line in the first facet, just below it on the right
        lab = ch[label == levels(d$label)[1]]
        p = p + geom_text(data = lab, aes(x = 3.45, y = chance - 0.03, label = 'random chance'),
                          inherit.aes = FALSE, hjust = 1, vjust = 1, size = 4, color = 'grey30', fontface = 'italic')
    }
    p
}
tum = s[, .(p_pick = mean(p_pick), n_events = .N), by = .(pair, assay_short)]
p_tum = ggplot(tum, aes(assay_short, p_pick, color = assay_short)) +
    geom_sina(size = 0.7, alpha = 0.6, stroke = 0) +
    geom_boxplot(aes(fill = assay_short), outlier.shape = NA, alpha = 0.12, width = 0.55, linewidth = 0.45) +
    scale_color_manual(values = assay_cols) + scale_fill_manual(values = assay_cols) +
    scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = expansion(add = 0.02)) +
    labs(x = NULL, y = 'Sample avg. identifiability') + theme_fig4
## D: fraction of events per type with P(MLE = truth) >= ident_cutoff; dotted line splits simple / complex
frac = s[type %in% c(simple_types, complex_types), .(frac = mean(p_pick >= ident_cutoff), n = .N), by = .(type, assay_short)]
frac[, type := factor(type, rev(c(simple_types, complex_types)))]
message('\nFraction of events with P(MLE = truth) >= ', ident_cutoff, ':')
print(dcast(frac, type ~ assay_short, value.var = 'frac')[order(-type)], digits = 2)
p_frac = ggplot(frac, aes(frac, type, fill = assay_short)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.8) +
    geom_hline(yintercept = length(complex_types) + 0.5, linetype = 'dotted', color = 'grey40', linewidth = 0.4) +
    scale_fill_manual(values = assay_cols, name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(add = c(0, 0.02))) +
    labs(x = bquote('Frac. events identifiable (P' >= .(ident_cutoff) * ')'), y = NULL, title = 'Identifiable events') +
    theme_fig4 + theme(panel.grid.major.y = element_blank(),
                       panel.grid.major.x = element_line(color = 'grey90', linewidth = 0.3),
                       legend.position = 'inside', legend.position.inside = c(0.98, 0.02),
                       legend.justification = c(1, 0), legend.key.size = unit(0.8, 'lines'),
                       plot.title = element_text(size = 13, face = 'bold', hjust = 0.5, color = 'grey10'))
## egg aligns the panel areas: A|C and B|D share rows, A/B and C/D share columns.
## widths/heights are panel sizes: type facets are 0.75 as wide as tall, C and D square
n_facets = max(length(simple_types), length(complex_types))
fig = egg::ggarrange(type_panel(simple_types, chance_label = TRUE), p_tum,
                     type_panel(complex_types), p_frac,
                     ncol = 2, widths = c(n_facets * 0.75, 1), heights = c(1, 1),
                     labels = c('A', 'C', 'B', 'D'), label.args = list(gp = grid::gpar(fontface = 'bold', fontsize = 14)),
                     draw = FALSE)
save_plot(fig, 'identifiability_figure', width = 13.5, height = 7)
