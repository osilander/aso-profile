# Fig S3 — Unfiltered vs specificity-filtered pooled coverage per gene
# pooled all-population greedy coverage
# Output: submission/figS2_coverage_filtered_vs_all.pdf

library(ggplot2)

K_MAX <- 10

# Unfiltered: same source as fig4A pooled black line
unfilt <- read.csv("gene_plot_inputs/coverage_curves_af0.1.csv",
                   stringsAsFactors = FALSE)
unfilt <- unfilt[unfilt$k <= K_MAX, c("k","frac_covered","gene")]
unfilt$filter <- "All common SNPs"

# Filtered: specificity-filtered pooled greedy
filt <- read.csv("figures/coverage_specificity.csv", stringsAsFactors = FALSE)
filt <- filt[filt$filter == "Pass specificity filters" & filt$k <= K_MAX,
             c("k","frac_covered","gene","filter")]

spec <- rbind(unfilt, filt)

# Gene order: match fig3a / fig4a
gene_ord <- c("CARD11","NLRP3","ABCC8","KIF1A","SCN2A","SCN8A","AKT1","KCNQ2",
              "STAT1","COL1A1","COL1A2","COL2A1","STAT3","SPTAN1","SPTLC1",
              "MTOR","KIF5B","KCNJ11","PIK3CA","RELA")
extra    <- setdiff(unique(spec$gene), gene_ord)
gene_ord <- c(gene_ord, extra)
spec$gene   <- factor(spec$gene,   levels = gene_ord)
spec$filter <- factor(spec$filter,
                      levels = c("All common SNPs", "Pass specificity filters"))

LINE_COLS <- c("All common SNPs"          = "grey55",
               "Pass specificity filters"  = "#0072B2")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor    = element_blank(),
        panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_line(colour = "grey88", linewidth = 0.3),
        strip.background    = element_rect(fill = "grey94", colour = "grey70"),
        strip.text          = element_text(size = 7, face = "bold"),
        legend.position     = "bottom",
        legend.title        = element_blank(),
        legend.key.width    = unit(1, "cm"))

p <- ggplot(spec, aes(x = k, y = frac_covered * 100,
                      colour = filter, group = filter)) +
  geom_line(linewidth = 0.7) +
  geom_point(shape = 21, fill = "white", size = 0.9, stroke = 0.5) +
  facet_wrap(~ gene, ncol = 5) +
  scale_colour_manual(values = LINE_COLS,
                      labels = c("All common SNPs"          = "All common SNPs         ",
                                 "Pass specificity filters" = "Common SNPs with no exact off-target binding site")) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  scale_y_continuous(limits = c(35, 101), breaks = seq(40, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Common SNP ASO targets added",
       y = "Individuals with ≥1 heterozygous target") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
  th

dir.create("submission", showWarnings = FALSE)
cairo_pdf("submission/figS3.pdf", width = 10, height = 8)
print(p)
dev.off()
cat("Written: submission/figS3.pdf\n")
