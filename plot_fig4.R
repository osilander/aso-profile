# Fig 4 — Common SNP ASO coverage and off-target specificity
#   A: Pooled all-population greedy coverage curves (all common SNPs, before filtering)
#   B: Per-gene before/after coverage at k=10 (pooled, horizontal dotplot)
#   C: Why filtering removes candidates — % no exact off-target + % repeat-overlapping
# Output: submission/fig4.pdf

library(ggplot2)
library(cowplot)
library(ggrepel)

GRP_COLS <- c(
  "Common SNP (minor)" = "#0072B2",
  "Common SNP (major)" = "#56B4E9",
  "Pathogenic"         = "#CC79A7",
  "Pathogenic (WT)"    = "#009E73"
)
GRP_ORDER <- names(GRP_COLS)
K_MAX  <- 10
K_SHOW <- 5

th_base <- theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(),
        legend.key.size  = unit(0.4, "cm"))

# ============================================================ #
# Panel A: pooled greedy coverage curves — all common SNPs
# ============================================================ #

grdy <- read.csv("gene_plot_inputs/coverage_curves_af0.1.csv",
                 stringsAsFactors = FALSE)
grdy <- grdy[grdy$k <= K_MAX, ]

# Order genes by coverage at k=10 (or max available), best to worst
ord10 <- grdy[grdy$k == K_MAX, ]
gene_ord <- ord10$gene[order(-ord10$frac_covered)]
# include genes that may plateau before k=10
missing <- setdiff(unique(grdy$gene), gene_ord)
gene_ord <- c(gene_ord, missing)
grdy$gene <- factor(grdy$gene, levels = gene_ord)

SUPERPOP_COLS <- c(AFR = "#E69F00", AMR = "#56B4E9",
                   EAS = "#009E73", EUR = "#0072B2", SAS = "#CC79A7")

pop_grdy <- read.csv("gene_plot_inputs/pop_greedy_curves_af0.1.csv",
                     stringsAsFactors = FALSE)
pop_grdy <- pop_grdy[pop_grdy$k <= K_MAX, ]
pop_grdy$gene    <- factor(pop_grdy$gene,    levels = gene_ord)
pop_grdy$superpop <- factor(pop_grdy$superpop, levels = names(SUPERPOP_COLS))

th_facet <- th_base +
  theme(strip.background  = element_rect(fill = "grey94", colour = "grey70"),
        strip.text         = element_text(size = 8.5, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
        legend.position    = "bottom",
        legend.spacing.x   = unit(0.3, "cm"))

pA <- ggplot(pop_grdy,
             aes(x = k, y = frac_covered * 100,
                 colour = superpop, group = superpop)) +
  geom_line(linewidth = 0.65, alpha = 0.8) +
  geom_point(shape = 21, fill = "white", size = 0.85, stroke = 0.45) +
  geom_line(data = grdy,
            aes(x = k, y = frac_covered * 100, group = gene),
            colour = "black", linewidth = 0.35, inherit.aes = FALSE) +
  geom_point(data = grdy,
             aes(x = k, y = frac_covered * 100, group = gene),
             colour = "black", shape = 21, fill = "white",
             size = 0.6, stroke = 0.35, inherit.aes = FALSE) +
  facet_wrap(~ gene, ncol = 5) +
  scale_colour_manual(values = SUPERPOP_COLS, name = NULL,
                      labels = c(AFR = "AFR            ", AMR = "AMR            ",
                                 EAS = "EAS            ", EUR = "EUR            ",
                                 SAS = "SAS")) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  scale_y_continuous(limits = c(35, 101), breaks = seq(40, 100, 10),
                     labels = function(x) ifelse(x %% 20 == 0, paste0(x, "%"), "")) +
  labs(x = "Common SNP targets added (k)",
       y = "Individuals with ≥1 heterozygous target") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
  th_facet

# ============================================================ #
# Panel B: before/after at k=10 (horizontal dotplot)
# ============================================================ #

spec <- read.csv("figures/coverage_specificity.csv", stringsAsFactors = FALSE)

get_plateau <- function(df_filter) {
  do.call(rbind, lapply(split(df_filter, df_filter$gene), function(d) {
    k_use <- if (K_SHOW %in% d$k) K_SHOW else max(d$k)
    data.frame(gene = d$gene[1], frac = d$frac_covered[d$k == k_use][1])
  }))
}

before <- get_plateau(spec[spec$filter == "All common SNPs", ])
after  <- get_plateau(spec[spec$filter == "Pass specificity filters", ])
names(before)[2] <- "before"
names(after)[2]  <- "after"

slope <- merge(before, after, by = "gene", all.x = TRUE)
slope$after[is.na(slope$after)] <- 0

# Sort by filtered coverage descending (top = best)
slope <- slope[order(-slope$after), ]
slope$gene <- factor(slope$gene, levels = slope$gene)

x_levels <- c("All common SNPs", "No exact\noff-target")

long_b <- rbind(
  data.frame(gene = slope$gene, x = x_levels[1], pct = slope$before * 100),
  data.frame(gene = slope$gene, x = x_levels[2], pct = slope$after  * 100)
)
long_b$x <- factor(long_b$x, levels = x_levels)

pB <- ggplot() +
  geom_hline(yintercept = c(80, 90, 95), linetype = "dashed",
             colour = "grey85", linewidth = 0.3) +
  geom_line(data = long_b,
            aes(x = x, y = pct, group = gene),
            colour = "grey60", linewidth = 0.45) +
  geom_point(data = long_b[long_b$x == x_levels[1], ],
             aes(x = x, y = pct),
             shape = 1, colour = "grey45", size = 2.2) +
  geom_point(data = long_b[long_b$x == x_levels[2], ],
             aes(x = x, y = pct),
             shape = 16, colour = "#0072B2", size = 2.5) +
  geom_text_repel(data = long_b[long_b$x == x_levels[2], ],
                  aes(x = x, y = pct, label = gene),
                  direction    = "y", hjust = 0, nudge_x = 0.15,
                  size = 2.5, segment.size = 0.25,
                  segment.colour = "grey65", box.padding = 0.15,
                  point.padding = 0.1, colour = "grey20",
                  max.overlaps = Inf, seed = 42) +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.7))) +
  scale_y_continuous(limits = c(40, 100), breaks = seq(40, 100, 10),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = sprintf("Individuals with ≥1 heterozygous target (k = %d)", K_SHOW)) +
  th_base +
  theme(panel.grid.major.x = element_blank(),
        legend.position    = "none")

# ============================================================ #
# Panel C: off-target pass rate + repeat overlap by class
# ============================================================ #

d <- read.delim("results/offtarget/offtarget_summary.tsv",
                stringsAsFactors = FALSE)
d$repeat_overlap_any <- as.logical(d$repeat_overlap_any)

d$group <- factor(with(d, ifelse(
  variant_class == "common_snp"                 & allele == "ALT", "Common SNP (minor)",
  ifelse(variant_class == "common_snp"          & allele == "REF", "Common SNP (major)",
  ifelse(variant_class == "pathogenic_missense"  & allele == "ALT", "Pathogenic",
         "Pathogenic (WT)")))), levels = GRP_ORDER)

sub_c <- d[d$index == "premrna", ]

pass_df <- aggregate(I(n_offtarget_0mm == 0) ~ group, data = sub_c, FUN = mean)
rep_df  <- aggregate(repeat_overlap_any  ~ group, data = sub_c, FUN = mean)
names(pass_df)[2] <- "val"
names(rep_df)[2]  <- "val"
pass_df$metric <- "No exact off-target hit"
rep_df$metric  <- "Repeat-overlapping"

bar_df <- rbind(pass_df, rep_df)
bar_df$val   <- bar_df$val * 100
bar_df$group <- factor(bar_df$group, levels = GRP_ORDER)
bar_df$metric <- factor(bar_df$metric,
                        levels = c("Repeat-overlapping", "No exact off-target hit"))

th_bar <- th_base +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5),
        strip.background = element_rect(fill = "grey94", colour = "grey70"),
        strip.text = element_text(size = 8),
        legend.position = "none")

pC <- ggplot(bar_df, aes(x = group, y = val, fill = group)) +
  geom_col(width = 0.55, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f%%", val)),
            vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = GRP_COLS) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  facet_wrap(~ metric) +
  labs(x = NULL, y = "% of ASO candidates\n(pre-mRNA index)") +
  th_bar

# ============================================================ #
# Assemble
# ============================================================ #

pA <- pA + theme(plot.margin = margin(4, 4, 8, 4, "mm"))
pB <- pB + theme(plot.margin = margin(8, 4, 4, 4, "mm"))
pC <- pC + theme(plot.margin = margin(8, 4, 4, 4, "mm"))

bot <- plot_grid(pB, pC, nrow = 1, labels = c("b", "c"),
                 rel_widths = c(1, 1.4),
                 label_size = 17, label_fontface = "bold")

fig4 <- plot_grid(pA, bot, nrow = 2, labels = c("a", ""),
                  rel_heights = c(1.9, 1),
                  label_size = 17, label_fontface = "bold")

dir.create("submission", showWarnings = FALSE)
cairo_pdf("submission/fig4.pdf", width = 10, height = 12)
print(fig4)
dev.off()
cat("Written: submission/fig4.pdf\n")
