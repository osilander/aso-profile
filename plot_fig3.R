# Fig 3 — Allele selectivity of common SNP ASO targets vs pathogenic SNVs
# Panel A: per-gene selectivity ECDFs (best_ddg for common SNPs; aso_alt for path)
# Panel B: mean ΔΔG asymmetry by substitution type (design strategy rationale)
# Output: figures/aso_selectivity.pdf

library(ggplot2)
library(cowplot)

# --------------------------------------------------------------------------- #
# Shared palette and theme
# --------------------------------------------------------------------------- #

CLASS_COLS <- c("Common SNP" = "#0072B2", "Pathogenic SNV" = "#CC79A7")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor  = element_blank(),
        legend.position   = "none",
        plot.margin       = margin(10, 12, 10, 10, "pt"))

# --------------------------------------------------------------------------- #
# Load and prepare data
# --------------------------------------------------------------------------- #

dat <- read.delim("results/aso_17mer_scores_extended.tsv",
                  stringsAsFactors = FALSE, check.names = FALSE)

dat$class <- ifelse(dat$variant_class == "common_snp",
                    "Common SNP", "Pathogenic SNV")

ok  <- dat[dat$pass_basic == 1, ]
ok$best_ddg <- pmax(ok$duplex_ddg_aso_ref, ok$duplex_ddg_aso_alt, na.rm = TRUE)

snp <- ok[ok$variant_class == "common_snp", ]
plp <- ok[ok$variant_class == "pathogenic_missense", ]

# --------------------------------------------------------------------------- #
# Panel A: Per-gene selectivity ECDFs (small multiples) — 4 lines per panel
# minor:  common SNP ASO targeting ALT (minor) allele  → duplex_ddg_aso_alt
# major:  common SNP ASO targeting REF (major) allele  → duplex_ddg_aso_ref
# path:   pathogenic ASO targeting ALT (pathogenic)    → duplex_ddg_aso_alt
# ref/WT: pathogenic ASO targeting REF (WT) allele     → duplex_ddg_aso_ref
# Gene order: descending median of common SNP minor-allele selectivity
# --------------------------------------------------------------------------- #

LINE_COLS <- c(
  "Common SNP (minor)"    = "#0072B2",
  "Common SNP (major)"    = "#56B4E9",
  "Pathogenic"            = "#CC79A7",
  "Pathogenic (WT)"       = "#009E73"
)

ecdf_dat <- rbind(
  data.frame(gene  = snp$gene, ddg = snp$duplex_ddg_aso_alt,
             grp = "Common SNP (minor)", stringsAsFactors = FALSE),
  data.frame(gene  = snp$gene, ddg = snp$duplex_ddg_aso_ref,
             grp = "Common SNP (major)", stringsAsFactors = FALSE),
  data.frame(gene  = plp$gene, ddg = plp$duplex_ddg_aso_alt,
             grp = "Pathogenic",          stringsAsFactors = FALSE),
  data.frame(gene  = plp$gene, ddg = plp$duplex_ddg_aso_ref,
             grp = "Pathogenic (WT)",    stringsAsFactors = FALSE)
)
ecdf_dat <- ecdf_dat[!is.na(ecdf_dat$ddg), ]

# Gene order: match fig4a (descending pooled greedy coverage at k=10)
gene_ord <- c("CARD11","NLRP3","ABCC8","KIF1A","SCN2A","SCN8A","AKT1","KCNQ2",
              "STAT1","COL1A1","COL1A2","COL2A1","STAT3","SPTAN1","SPTLC1",
              "MTOR","KIF5B","KCNJ11","PIK3CA","RELA")
# keep any genes in data not in the list (future-proof)
extra <- setdiff(unique(ecdf_dat$gene), gene_ord)
gene_ord <- c(gene_ord, extra)
ecdf_dat$gene_f <- factor(ecdf_dat$gene, levels = gene_ord)
ecdf_dat$grp    <- factor(ecdf_dat$grp,  levels = names(LINE_COLS))

pA <- ggplot(ecdf_dat, aes(x = ddg, colour = grp, group = grp)) +
  stat_ecdf(linewidth = 0.55) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3,
             linetype = "dashed") +
  facet_wrap(~ gene_f, ncol = 5) +
  scale_colour_manual(values = LINE_COLS, name = NULL) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(x = "Selectivity ΔΔG (kcal/mol)  →  more selective", y = "Cumulative fraction",
       title = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background   = element_rect(fill = "grey94", colour = "grey70"),
        strip.text         = element_text(size = 7.5, face = "bold"),
        legend.position    = "none",
        plot.margin        = margin(10, 12, 10, 10, "pt"))

# --------------------------------------------------------------------------- #
# Panel B: Mean ΔΔG asymmetry by substitution type × variant class
# Shows design strategy: which allele to target for best selectivity
# --------------------------------------------------------------------------- #

ok$sub_type <- paste0(toupper(ok$ref), ">", toupper(ok$alt))
ok$delta    <- ok$duplex_ddg_aso_ref - ok$duplex_ddg_aso_alt

PAIRS <- list(
  "A–G" = c("A>G", "G>A"), "C–T" = c("C>T", "T>C"),
  "A–C" = c("A>C", "C>A"), "C–G" = c("C>G", "G>C"),
  "A–T" = c("A>T", "T>A"), "G–T" = c("G>T", "T>G")
)
pair_map <- do.call(rbind, lapply(names(PAIRS), function(p)
  data.frame(sub_type = PAIRS[[p]], pair = p, stringsAsFactors = FALSE)))

agg_d   <- aggregate(delta ~ sub_type + variant_class, data = ok, FUN = mean, na.rm = TRUE)
agg_sd  <- aggregate(delta ~ sub_type + variant_class, data = ok, FUN = sd,   na.rm = TRUE)
agg_n   <- aggregate(delta ~ sub_type + variant_class, data = ok, FUN = length)
colnames(agg_sd)[3] <- "sd"
colnames(agg_n)[3]  <- "n"
bar_dat <- merge(merge(merge(agg_d, agg_sd, by = c("sub_type", "variant_class")),
                       agg_n,    by = c("sub_type", "variant_class")),
                 pair_map, by = "sub_type")
bar_dat$pair         <- factor(bar_dat$pair, levels = names(PAIRS))
bar_dat$variant_class <- factor(bar_dat$variant_class,
  levels = c("common_snp", "pathogenic_missense"),
  labels = c("Common SNP", "Pathogenic SNV"))
bar_dat$sub_type  <- factor(bar_dat$sub_type)
bar_dat$hjust_n   <- ifelse(bar_dat$delta >= 0, -0.15, 1.15)

pB <- ggplot(bar_dat,
             aes(y = sub_type, x = delta, fill = variant_class)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, alpha = 0.9) +
  geom_text(aes(label = n,
                x     = delta + sign(delta) * 0.1,
                hjust = hjust_n),
            position = position_dodge(width = 0.72),
            size = 2.6, colour = "grey30") +
  scale_fill_manual(values = CLASS_COLS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(y        = NULL,
       x        = "Mean ΔΔG asymmetry (kcal/mol)",
       subtitle = "← ALT / path allele more selective                 REF / WT more selective →",
       title    = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position    = "none",
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.placement    = "outside",
        strip.text.y.left  = element_text(size = 8.5, face = "bold", angle = 0),
        strip.background   = element_rect(fill = "grey94", colour = "grey70"),
        plot.subtitle      = element_text(size = 6.5, colour = "grey45",
                                          face = "italic", hjust = 0.5,
                                          margin = margin(b = 4)),
        plot.margin        = margin(10, 12, 10, 10, "pt"))

# --------------------------------------------------------------------------- #
# Shared legend + assemble
# --------------------------------------------------------------------------- #

leg_dat <- data.frame(grp = factor(names(LINE_COLS), levels = names(LINE_COLS)),
                      x = 0, y = 1)
legend_src <- ggplot(leg_dat, aes(x = x, y = y, colour = grp, group = grp)) +
  geom_line() +
  scale_colour_manual(values = LINE_COLS, name = NULL,
                      labels = c("Common SNP (minor)" = "Common SNP (minor)         ",
                                 "Common SNP (major)" = "Common SNP (major)",
                                 "Pathogenic"         = "Pathogenic                 ",
                                 "Pathogenic (WT)"    = "Pathogenic (WT)")) +
  guides(colour = guide_legend(nrow = 2,
                               override.aes = list(linewidth = 1.4))) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 9),
        legend.key.width = unit(0.5, "cm"))
shared_legend <- get_legend(legend_src)

# separate fill legend for panel B
legend_B <- get_legend(
  ggplot(bar_dat, aes(y = sub_type, x = delta, fill = variant_class)) +
    geom_col() +
    scale_fill_manual(values = CLASS_COLS, name = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 9),
          legend.key.size = unit(0.45, "cm"))
)

panels <- plot_grid(pA, pB,
                    ncol           = 2,
                    rel_widths     = c(2, 1),
                    labels         = c("a", "b"),
                    label_size     = 14,
                    label_fontface = "bold",
                    align          = "hv",
                    axis           = "tblr")

legends <- plot_grid(shared_legend, legend_B, ncol = 2, rel_widths = c(2, 1))

dir.create("figures", showWarnings = FALSE)
cairo_pdf("figures/aso_selectivity.pdf", width = 10, height = 5.5, onefile = FALSE)
plot_grid(panels, legends, nrow = 2, rel_heights = c(1, 0.1))
dev.off()
cat("Written: figures/aso_selectivity.pdf\n")
