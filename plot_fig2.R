# Fig 2 — ASO design context: RNA accessibility and 4-mer score
# Panel A: % unpaired at 17-mer (RNAfold)
# Panel B: 17-mer sequence score (mean screen-adjusted 4-mer coef), per gene
# Output: figures/aso_context.pdf

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

dat <- read.delim("aso_17mer_scores_extended.tsv",
                  stringsAsFactors = FALSE, check.names = FALSE)

dat$class <- ifelse(dat$variant_class == "common_snp",
                    "Common SNP", "Pathogenic SNV")

ok <- dat[dat$pass_basic == 1, ]

# --------------------------------------------------------------------------- #
# 17-mer sequence score: mean screen-adjusted 4-mer coefficient
# --------------------------------------------------------------------------- #

fourmer <- read.csv("data/all_4mer_motif_effects_augmented_aso_atlas.csv",
                    stringsAsFactors = FALSE)
coef_map <- setNames(fourmer$screen_adjusted_coef, fourmer$motif)

score_17mer <- function(seq) {
  kmers <- sapply(1:(nchar(seq) - 3), function(i) substr(seq, i, i + 3))
  mean(coef_map[kmers], na.rm = TRUE)
}
ok$score_17mer <- sapply(ok$seq_17mer, score_17mer)

# --------------------------------------------------------------------------- #
# Parse RNAfold dot-bracket output → fraction unpaired at 17-mer window
# --------------------------------------------------------------------------- #

parse_unpaired <- function(rnafold_file) {
  lines <- readLines(rnafold_file)
  h_idx <- which(startsWith(lines, ">"))
  do.call(rbind, lapply(h_idx, function(i) {
    parts    <- strsplit(sub("^>", "", lines[i]), "\\|")[[1]]
    join_key <- paste(parts[1:3], collapse = "|")
    struct   <- sub(" \\(.*\\)$", "", lines[i + 2])
    w17      <- substr(struct, 93, 109)
    data.frame(join_key = join_key,
               frac_unpaired_17mer = mean(strsplit(w17, "")[[1]] == "."),
               stringsAsFactors = FALSE)
  }))
}

cat("Parsing RNAfold outputs ...\n")
unp <- rbind(
  parse_unpaired("results/rnafold/common_snps.plusminus100.rnafold.txt"),
  parse_unpaired("results/rnafold/pathogenic_missense.plusminus100.rnafold.txt")
)

ok$join_key <- paste(ok$variant_class, ok$gene, ok$variant_id, sep = "|")
ok <- merge(ok, unp, by = "join_key", all.x = TRUE)

# --------------------------------------------------------------------------- #
# Panel A: Fraction of 17-mer bases unpaired
# --------------------------------------------------------------------------- #

ok_A <- ok[!is.na(ok$frac_unpaired_17mer), ]
ok_A$wt <- 1 / tabulate(factor(ok_A$class))[factor(ok_A$class)]

pA <- ggplot(ok_A, aes(x = frac_unpaired_17mer * 100, weight = wt,
                       fill = class, colour = class)) +
  geom_histogram(position = "identity", binwidth = 100 / 17, alpha = 0.55,
                 linewidth = 0.3) +
  scale_fill_manual(values = CLASS_COLS, name = NULL,
                    labels = c("Common SNP" = "Common SNP        ",
                               "Pathogenic SNV" = "Pathogenic SNV")) +
  scale_colour_manual(values = CLASS_COLS, name = NULL,
                      labels = c("Common SNP" = "Common SNP        ",
                                 "Pathogenic SNV" = "Pathogenic SNV")) +
  scale_x_continuous(limits = c(-3, 103),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Target site unpaired (%)", y = "Fraction", title = NULL) +
  annotate("text", x = 100, y = Inf, label = "more accessible to ASO binding →",
           hjust = 1, vjust = 1.8, size = 2.7, colour = "black",
           fontface = "italic") +
  th

# --------------------------------------------------------------------------- #
# Panel B: 17-mer sequence score — overlaid histogram
# --------------------------------------------------------------------------- #

ok_B <- ok[!is.na(ok$score_17mer), ]
ok_B$wt <- 1 / tabulate(factor(ok_B$class))[factor(ok_B$class)]

pB <- ggplot(ok_B, aes(x = score_17mer, weight = wt,
                       fill = class, colour = class)) +
  geom_histogram(position = "identity", binwidth = 0.3125, alpha = 0.55,
                 linewidth = 0.3) +
  scale_fill_manual(values = CLASS_COLS, name = NULL,
                    labels = c("Common SNP" = "Common SNP        ",
                               "Pathogenic SNV" = "Pathogenic SNV")) +
  scale_colour_manual(values = CLASS_COLS, name = NULL,
                      labels = c("Common SNP" = "Common SNP        ",
                                 "Pathogenic SNV" = "Pathogenic SNV")) +
  labs(x = "4-mer sequence score", y = "Fraction", title = NULL) +
  annotate("text", x = -Inf, y = Inf, label = "increased mRNA knockdown efficacy →",
           hjust = -0.05, vjust = 1.4, size = 2.7, colour = "black",
           fontface = "italic") +
  th

# --------------------------------------------------------------------------- #
# Shared legend + assemble
# --------------------------------------------------------------------------- #

legend_src <- ggplot(data.frame(class = factor(names(CLASS_COLS),
                                               levels = names(CLASS_COLS)),
                                x = 0, y = 1),
                     aes(x = x, y = y, fill = class)) +
  geom_col() +
  scale_fill_manual(values = CLASS_COLS, name = NULL,
                    labels = c("Common SNP" = "Common SNP        ",
                               "Pathogenic SNV" = "Pathogenic SNV")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))
shared_legend <- get_legend(legend_src)

dir.create("figures", showWarnings = FALSE)
cairo_pdf("figures/aso_context.pdf", width = 5.04, height = 5.6, onefile = FALSE)
plot_grid(pA, shared_legend, pB,
          nrow           = 3,
          labels         = c("a", "", "b"),
          label_size     = 10,
          label_fontface = "bold",
          rel_heights    = c(1, 0.12, 1),
          align          = "v",
          axis           = "lr")
dev.off()
cat("Written: figures/aso_context.pdf\n")
