# Fig S2 (ancestry-stratified filtered) + Fig S4 (random vs greedy per-pop)
# Output: figures/ancestry_coverage_af0.1.pdf
#         figures/ancestry_coverage_specificity_af0.1.pdf

library(ggplot2)
library(cowplot)

SUPERPOP_COLS <- c(AFR = "#E69F00", AMR = "#56B4E9", EAS = "#009E73",
                   EUR = "#0072B2", SAS = "#CC79A7")
FILTER_COLS <- c(
  "All common SNPs"          = "grey70",
  "Pass basic filters"       = "#56B4E9",
  "Pass specificity filters" = "#CC79A7"
)
K_MAX <- 10

th_facet <- theme_bw(base_size = 9) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background   = element_rect(fill = "grey94", colour = "grey70"),
        strip.text         = element_text(size = 8, face = "bold"),
        legend.position    = "bottom",
        legend.key.width   = unit(1.2, "cm"),
        legend.spacing.x   = unit(0.4, "cm"))

# ============================================================ #
# VERSION 1: all common SNPs, ancestry-stratified
# ============================================================ #

anc_df      <- read.csv("gene_plot_inputs/ancestry_coverage_curves_af0.1.csv",
                        stringsAsFactors = FALSE)
grdy_df     <- read.csv("gene_plot_inputs/pop_greedy_curves_af0.1.csv",
                        stringsAsFactors = FALSE)
rand_all_df <- read.csv("gene_plot_inputs/random_all_curves_af0.1.csv",
                        stringsAsFactors = FALSE)
greedy_df   <- read.csv("gene_plot_inputs/coverage_curves_af0.1.csv",
                        stringsAsFactors = FALSE)

anc_df      <- anc_df[anc_df$k      <= K_MAX, ]
grdy_df     <- grdy_df[grdy_df$k    <= K_MAX, ]
rand_all_df <- rand_all_df[rand_all_df$k <= K_MAX, ]
greedy_df   <- greedy_df[greedy_df$k    <= K_MAX, ]

gene_levels <- c("CARD11","NLRP3","ABCC8","KIF1A","SCN2A","SCN8A","AKT1","KCNQ2",
                 "STAT1","COL1A1","COL1A2","COL2A1","STAT3","SPTAN1","SPTLC1",
                 "MTOR","KIF5B","KCNJ11","PIK3CA","RELA")

anc_df$gene      <- factor(anc_df$gene,     levels = gene_levels)
anc_df$superpop  <- factor(anc_df$superpop, levels = names(SUPERPOP_COLS))
grdy_df$gene     <- factor(grdy_df$gene,    levels = gene_levels)
grdy_df$superpop <- factor(grdy_df$superpop, levels = names(SUPERPOP_COLS))
greedy_df$gene   <- factor(greedy_df$gene,  levels = gene_levels)

p1 <- ggplot(anc_df,
             aes(x = k, y = mean * 100, colour = superpop, group = superpop)) +
  geom_hline(yintercept = 95, linetype = "dashed",
             colour = "grey70", linewidth = 0.35) +
  geom_line(linewidth = 0.5, alpha = 0.45, linetype = "dashed") +
  geom_line(data = grdy_df,
            aes(x = k, y = frac_covered * 100, colour = superpop, group = superpop),
            linewidth = 0.9, linetype = "solid", alpha = 0.7,
            inherit.aes = FALSE) +
  geom_point(data = grdy_df,
             aes(x = k, y = frac_covered * 100, colour = superpop, group = superpop),
             shape = 21, fill = "white", size = 1.0, stroke = 0.55,
             inherit.aes = FALSE) +
  geom_line(data = greedy_df,
            aes(x = k, y = frac_covered * 100, group = gene),
            colour = "black", linewidth = 0.4,
            inherit.aes = FALSE) +
  geom_point(data = greedy_df,
             aes(x = k, y = frac_covered * 100, group = gene),
             colour = "black", shape = 21, fill = "white",
             size = 0.7, stroke = 0.45,
             inherit.aes = FALSE) +
  facet_wrap(~ gene, ncol = 5) +
  scale_colour_manual(values = SUPERPOP_COLS, name = "Optimal (by population)",
                      labels = c(AFR = "AFR            ", AMR = "AMR            ",
                                 EAS = "EAS            ", EUR = "EUR            ",
                                 SAS = "SAS")) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  scale_y_continuous(limits = c(20, 101),
                     breaks = seq(20, 100, by = 20),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Common SNP ASO targets added",
       y = "Individuals with ≥1 heterozygous ASO target") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8), order = 1)) +
  th_facet

dir.create("figures", showWarnings = FALSE)
cairo_pdf("figures/ancestry_coverage_af0.1.pdf", width = 10, height = 8)
print(p1)
dev.off()
cat("Written: figures/ancestry_coverage_af0.1.pdf\n")

# ============================================================ #
# VERSION 2: specificity-filtered SNPs, ancestry-stratified
# ============================================================ #

pop_spec_file <- "figures/pop_greedy_specificity.csv"
spec_df       <- read.csv("figures/coverage_specificity.csv", stringsAsFactors = FALSE)

if (file.exists(pop_spec_file)) {
  pop_spec_df <- read.csv(pop_spec_file, stringsAsFactors = FALSE)
  pop_spec_df <- pop_spec_df[pop_spec_df$k <= K_MAX, ]
  pop_spec_df$gene     <- factor(pop_spec_df$gene,     levels = gene_levels)
  pop_spec_df$superpop <- factor(pop_spec_df$superpop, levels = names(SUPERPOP_COLS))

  # Pooled specificity line (greedy, all individuals)
  pool_spec <- spec_df[spec_df$filter == "Pass specificity filters" &
                         spec_df$k <= K_MAX, ]
  pool_spec$gene <- factor(pool_spec$gene, levels = gene_levels)

  p2 <- ggplot(pop_spec_df,
               aes(x = k, y = frac_covered * 100,
                   colour = superpop, group = superpop)) +
    geom_hline(yintercept = 95, linetype = "dashed",
               colour = "grey70", linewidth = 0.35) +
    geom_line(linewidth = 0.9, linetype = "solid", alpha = 0.7) +
    geom_point(shape = 21, fill = "white", size = 1.0, stroke = 0.55) +
    geom_line(data = pool_spec,
              aes(x = k, y = frac_covered * 100, group = gene),
              colour = "black", linewidth = 0.4,
              inherit.aes = FALSE) +
    geom_point(data = pool_spec,
               aes(x = k, y = frac_covered * 100, group = gene),
               colour = "black", shape = 21, fill = "white",
               size = 0.7, stroke = 0.45,
               inherit.aes = FALSE) +
    facet_wrap(~ gene, ncol = 5) +
    scale_colour_manual(values = SUPERPOP_COLS, name = "Optimal (by population)",
                        labels = c(AFR = "AFR            ", AMR = "AMR            ",
                                   EAS = "EAS            ", EUR = "EUR            ",
                                   SAS = "SAS")) +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    scale_y_continuous(limits = c(0, 101),
                       breaks = seq(0, 100, by = 20),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Specificity-filtered ASO targets added",
         y = "Individuals with ≥1 heterozygous ASO target") +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.8), order = 1)) +
    th_facet

  cairo_pdf("figures/ancestry_coverage_specificity_af0.1.pdf", width = 10, height = 8)
  print(p2)
  dev.off()
  cat("Written: figures/ancestry_coverage_specificity_af0.1.pdf\n")

} else {
  # Fallback: three filter tier lines if pop-specific data not yet available
  warning("figures/pop_greedy_specificity.csv not found — run coverage_curves.R first")

  spec_df <- spec_df[spec_df$k <= K_MAX, ]
  spec_df$filter <- factor(spec_df$filter,
                           levels = c("All common SNPs",
                                      "Pass basic filters",
                                      "Pass specificity filters"))
  spec_df$gene   <- factor(spec_df$gene, levels = gene_levels)

  p2 <- ggplot(spec_df,
               aes(x = k, y = frac_covered * 100,
                   colour = filter, group = interaction(gene, filter))) +
    geom_hline(yintercept = 95, linetype = "dashed",
               colour = "grey70", linewidth = 0.35) +
    geom_line(linewidth = 0.75) +
    geom_point(shape = 21, fill = "white", size = 1.0, stroke = 0.55) +
    facet_wrap(~ gene, ncol = 5) +
    scale_colour_manual(values = FILTER_COLS, name = NULL) +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    scale_y_continuous(limits = c(0, 101),
                       breaks = seq(0, 100, by = 20),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Common SNP ASO targets added",
         y = "Individuals with ≥1 heterozygous ASO target") +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.8))) +
    th_facet

  cairo_pdf("figures/ancestry_coverage_specificity_af0.1.pdf", width = 10, height = 8)
  print(p2)
  dev.off()
  cat("Written: figures/ancestry_coverage_specificity_af0.1.pdf\n")
}
