# Tables S1 and S2 — coverage thresholds and candidate counts
# Table S1: k at 80/90/95% coverage before and after specificity filtering, per gene
# Table S2: candidate counts retained after each filter stage
# Output: submission/tableS1_coverage_thresholds.tsv
#         submission/tableS2_candidate_counts.tsv

dir.create("submission", showWarnings = FALSE)

# ============================================================ #
# Table S1: coverage thresholds
# ============================================================ #

spec <- read.csv("figures/coverage_specificity.csv", stringsAsFactors = FALSE)

k_thresh <- function(df, filter_name, threshold) {
  sub <- df[df$filter == filter_name & df$frac_covered >= threshold, ]
  if (nrow(sub) == 0) return(NA_integer_)
  min(sub$k)
}

genes <- sort(unique(spec$gene))
filters <- c("All common SNPs", "Pass specificity filters")
thresholds <- c(0.80, 0.90, 0.95)
thresh_labs <- c("k80", "k90", "k95")

rows <- lapply(genes, function(g) {
  d <- spec[spec$gene == g, ]
  r <- data.frame(gene = g)
  for (fi in seq_along(filters)) {
    prefix <- if (fi == 1) "unfiltered" else "filtered"
    for (ti in seq_along(thresholds)) {
      col <- paste0(prefix, "_", thresh_labs[ti])
      r[[col]] <- k_thresh(d, filters[fi], thresholds[ti])
    }
  }
  r
})
t1 <- do.call(rbind, rows)

# Order by unfiltered k95 then gene name
t1 <- t1[order(t1$unfiltered_k95, t1$gene, na.last = TRUE), ]

write.table(t1, "submission/tableS1_coverage_thresholds.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat("Written: submission/tableS1_coverage_thresholds.tsv\n")

# ============================================================ #
# Table S2: candidate counts per filter stage
# ============================================================ #

dat <- read.delim("results/aso_17mer_scores_extended.tsv",
                  stringsAsFactors = FALSE, check.names = FALSE)
ot  <- read.delim("results/offtarget/offtarget_summary.tsv",
                  stringsAsFactors = FALSE)
ot$pass_specificity_moderate <- as.logical(ot$pass_specificity_moderate)
ot$repeat_overlap_any        <- as.logical(ot$repeat_overlap_any)

# common SNPs: variant-level counts from main TSV
cs       <- dat[dat$variant_class == "common_snp", ]
cs_basic <- cs[cs$pass_basic == 1, ]

# Stage 1: all common SNP variants
n_all   <- length(unique(cs$variant_id))
# Stage 2: pass basic filters
n_basic <- length(unique(cs_basic$variant_id))

# Stage 3: pass specificity (moderate) — AND across indices, OR across alleles
ot_cs <- ot[ot$variant_class == "common_snp", ]

# per variant_id × allele: all across indices
ot_per_allele <- aggregate(pass_specificity_moderate ~ variant_id + allele,
                           data = ot_cs, FUN = all)
# per variant_id: any allele passes
ot_per_var <- aggregate(pass_specificity_moderate ~ variant_id,
                        data = ot_per_allele, FUN = any)
pass_var_ids <- ot_per_var$variant_id[ot_per_var$pass_specificity_moderate]

n_spec_allele  <- sum(ot_per_allele$pass_specificity_moderate)
n_spec_variant <- length(pass_var_ids)

# Per-gene breakdown using variant_id gene lookup from main TSV
var_gene <- unique(cs[, c("variant_id","gene")])

per_gene <- do.call(rbind, lapply(sort(unique(cs$gene)), function(g) {
  g_vids  <- cs$variant_id[cs$gene == g]
  g_all   <- length(unique(g_vids))
  g_basic <- length(unique(cs_basic$variant_id[cs_basic$gene == g]))

  # allele-level pass for this gene
  g_allele <- ot_per_allele[ot_per_allele$variant_id %in% g_vids, ]
  g_spec_a <- sum(g_allele$pass_specificity_moderate)
  g_spec_v <- sum(ot_per_var$pass_specificity_moderate[
                    ot_per_var$variant_id %in% g_vids])

  data.frame(gene = g, all_variants = g_all, pass_basic = g_basic,
             pass_specificity_alleles = g_spec_a,
             pass_specificity_variants = g_spec_v)
}))

totals <- data.frame(gene = "TOTAL",
                     all_variants = n_all, pass_basic = n_basic,
                     pass_specificity_alleles = n_spec_allele,
                     pass_specificity_variants = n_spec_variant)
t2 <- rbind(per_gene, totals)

write.table(t2, "submission/tableS2_candidate_counts.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Written: submission/tableS2_candidate_counts.tsv\n")
