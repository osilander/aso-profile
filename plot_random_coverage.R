#!/usr/bin/env Rscript
# Random-selection coverage curves (faceted by gene) vs greedy ceiling.
# Requires: vcf_cache/*_het_matrix.rds (sync from HPC vcf_cache/)
#           gene_plot_inputs/coverage_curves.csv (greedy curves)

library(ggplot2)

RDS_DIR  <- "vcf_cache"
N_SIM    <- 500
K_MAX    <- 10
SET_SEED <- 42

# --------------------------------------------------------------------------- #
# Simulate random coverage for one gene
# Returns data.frame: k, mean, lo (10th pct), hi (90th pct)
# --------------------------------------------------------------------------- #

sim_random <- function(het_mat, k_max = K_MAX, n_sim = N_SIM) {
  n_indiv <- ncol(het_mat)
  n_snps  <- nrow(het_mat)
  k_max   <- min(k_max, n_snps)

  mat <- matrix(NA_real_, n_sim, k_max)
  for (i in seq_len(n_sim)) {
    perm    <- sample(n_snps)
    covered <- rep(FALSE, n_indiv)
    for (k in seq_len(k_max)) {
      covered    <- covered | het_mat[perm[k], ]
      mat[i, k]  <- sum(covered) / n_indiv
    }
  }

  data.frame(
    k    = seq_len(k_max),
    mean = colMeans(mat),
    lo   = apply(mat, 2, quantile, 0.10),
    hi   = apply(mat, 2, quantile, 0.90)
  )
}

# --------------------------------------------------------------------------- #
# Run simulations across all cached genes
# --------------------------------------------------------------------------- #

set.seed(SET_SEED)
rds_files <- list.files(RDS_DIR, pattern = "_het_matrix\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No .rds files found in ", RDS_DIR)

rand_list <- lapply(rds_files, function(f) {
  gene <- sub("_het_matrix\\.rds$", "", basename(f))
  cat(sprintf("%-10s simulating ...\n", gene))
  mat  <- readRDS(f)
  d    <- sim_random(mat)
  d$gene <- gene
  d
})
rand_df <- do.call(rbind, rand_list)

# --------------------------------------------------------------------------- #
# Load greedy curves (already computed)
# --------------------------------------------------------------------------- #

greedy_df <- read.csv("gene_plot_inputs/coverage_curves.csv",
                      stringsAsFactors = FALSE)
greedy_df <- greedy_df[greedy_df$k <= K_MAX, ]

# Keep only genes present in both
genes <- intersect(unique(rand_df$gene), unique(greedy_df$gene))
rand_df   <- rand_df[rand_df$gene     %in% genes, ]
greedy_df <- greedy_df[greedy_df$gene %in% genes, ]

# Consistent gene ordering: by greedy coverage at k=5 descending
ord <- greedy_df[greedy_df$k == 5, ]
ord <- ord[order(-ord$frac_covered), "gene"]
gene_levels <- ord
rand_df$gene   <- factor(rand_df$gene,   levels = gene_levels)
greedy_df$gene <- factor(greedy_df$gene, levels = gene_levels)

# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #

p <- ggplot() +
  # 95% reference line
  geom_hline(yintercept = 95, linetype = "dashed",
             colour = "grey65", linewidth = 0.35) +
  # Random: shaded band (10th–90th pct)
  geom_ribbon(data = rand_df,
              aes(x = k, ymin = lo * 100, ymax = hi * 100),
              fill = "#4e79a7", alpha = 0.18) +
  # Random: mean line
  geom_line(data = rand_df,
            aes(x = k, y = mean * 100),
            colour = "#4e79a7", linewidth = 0.7, linetype = "dashed") +
  # Greedy: solid ceiling
  geom_line(data = greedy_df,
            aes(x = k, y = frac_covered * 100),
            colour = "#e15759", linewidth = 0.85) +
  geom_point(data = greedy_df,
             aes(x = k, y = frac_covered * 100),
             colour = "#e15759", size = 1.4) +
  facet_wrap(~ gene, ncol = 5) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  scale_y_continuous(limits = c(NA, 101),
                     breaks = c(50, 75, 95),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Common SNP handles added",
       y = "Individuals with ≥1 heterozygous handle") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background   = element_rect(fill = "grey94", colour = "grey70"),
        strip.text         = element_text(size = 8, face = "bold"),
        plot.margin        = margin(8, 8, 8, 8, "pt"))

dir.create("figures", showWarnings = FALSE)
cairo_pdf("figures/random_coverage.pdf", width = 10, height = 8)
print(p)
dev.off()
cat("Written: figures/random_coverage.pdf\n")
