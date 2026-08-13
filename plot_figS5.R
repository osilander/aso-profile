# Fig S5 — Pre-mRNA off-target summary (ECDFs + bar summaries)
# Input: results/offtarget/offtarget_summary.tsv
# Output: submission/figS3.pdf

library(ggplot2)
library(cowplot)

GRP_COLS <- c(
  "Common SNP (minor)"  = "#0072B2",
  "Common SNP (major)"  = "#56B4E9",
  "Pathogenic"          = "#CC79A7",
  "Pathogenic (WT)"     = "#009E73"
)
GRP_ORDER <- names(GRP_COLS)

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "none")

d <- read.delim("results/offtarget/offtarget_summary.tsv",
                stringsAsFactors = FALSE)

# Python writes True/False as strings; coerce to logical
d$repeat_overlap_any          <- as.logical(d$repeat_overlap_any)
d$pass_specificity_strict     <- as.logical(d$pass_specificity_strict)
d$pass_specificity_moderate   <- as.logical(d$pass_specificity_moderate)
d$pass_specificity_permissive <- as.logical(d$pass_specificity_permissive)

d$group <- factor(with(d, ifelse(
  variant_class == "common_snp"                 & allele == "ALT", "Common SNP (minor)",
  ifelse(variant_class == "common_snp"          & allele == "REF", "Common SNP (major)",
  ifelse(variant_class == "pathogenic_missense"  & allele == "ALT", "Pathogenic",
         "Pathogenic (WT)")))), levels = GRP_ORDER)

# premRNA only
d <- d[d$index == "premrna", ]

sub <- d

log_breaks <- c(1, 10, 100, 1000)

# ------------------------------------------------------------------ #
# A: ECDF of 0mm (exact) off-target count
# ------------------------------------------------------------------ #
pA <- ggplot(sub, aes(x = n_offtarget_0mm + 1,
                      colour = group, group = group)) +
  stat_ecdf(linewidth = 0.7) +
  scale_colour_manual(values = GRP_COLS) +
  scale_x_log10(breaks = log_breaks,
                labels = as.character(log_breaks)) +
  labs(x = "Exact off-target hits + 1",
       y = "Fraction of ASO targets",
       title = "Exact off-target hits") +
  th

# ------------------------------------------------------------------ #
# B: ECDF of 1mm off-target count
# ------------------------------------------------------------------ #
pB <- ggplot(sub, aes(x = n_offtarget_1mm + 1,
                      colour = group, group = group)) +
  stat_ecdf(linewidth = 0.7) +
  scale_colour_manual(values = GRP_COLS) +
  scale_x_log10(breaks = log_breaks,
                labels = as.character(log_breaks)) +
  labs(x = "1-mismatch off-target hits + 1",
       y = "Fraction of ASO targets",
       title = "1-mismatch off-target hits") +
  th

# ------------------------------------------------------------------ #
# C: % with no exact off-target hit (bar)
# ------------------------------------------------------------------ #
sub$no_exact <- sub$n_offtarget_0mm == 0
exact_df <- aggregate(no_exact ~ group, data = sub, FUN = mean)
exact_df$group <- factor(exact_df$group, levels = GRP_ORDER)

pC <- ggplot(exact_df, aes(x = group, y = no_exact * 100, fill = group)) +
  geom_col(width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = GRP_COLS) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "% with no exact off-target hit") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

# ------------------------------------------------------------------ #
# D: repeat overlap % (bar)
# ------------------------------------------------------------------ #
rep_df <- aggregate(repeat_overlap_any ~ group, data = sub, FUN = mean)
rep_df$group <- factor(rep_df$group, levels = GRP_ORDER)

pD <- ggplot(rep_df, aes(x = group, y = repeat_overlap_any * 100, fill = group)) +
  geom_col(width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = GRP_COLS) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "% target site overlapping repeat") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

# ------------------------------------------------------------------ #
# Legend
# ------------------------------------------------------------------ #
leg_df  <- data.frame(group = factor(GRP_ORDER, levels = GRP_ORDER), y = 1)
leg_src <- ggplot(leg_df, aes(x = group, y = y, colour = group)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = GRP_COLS, name = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 8),
        legend.key.width = unit(0.6, "cm"))
shared_legend <- get_legend(leg_src)

# ------------------------------------------------------------------ #
# Assemble: ECDFs top row (A, B), bars bottom row (C, D)
# ------------------------------------------------------------------ #
top <- plot_grid(pA, pB, nrow = 1,
                 labels = c("a", "b"),
                 label_size = 11, label_fontface = "bold")

bot <- plot_grid(pC, pD, nrow = 1,
                 labels = c("c", "d"),
                 label_size = 11, label_fontface = "bold")

panels <- plot_grid(top, bot, nrow = 2, rel_heights = c(1, 0.9))

dir.create("figures",        showWarnings = FALSE)
dir.create("submission", showWarnings = FALSE)
for (out in c("figures/offtarget_allele_groups_premrna.pdf",
              "submission/figS5.pdf")) {
  cairo_pdf(out, width = 8.1, height = 6.3)
  print(plot_grid(panels, shared_legend, nrow = 2, rel_heights = c(1, 0.07)))
  dev.off()
  cat("Written:", out, "\n")
}
