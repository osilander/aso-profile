#!/usr/bin/env Rscript
# Greedy set-cover coverage curves: for each gene, how many optimally chosen
# common SNP handles are needed so that X% of 1KGP individuals are heterozygous
# for at least one?
#
# Data source: 1KGP Phase 3, GRCh38, phased biallelic SNVs
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/
#   1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

library(ggplot2)
library(dplyr)

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

TSV_FILE      <- "results/aso_17mer_scores_extended.tsv"
CACHE_DIR     <- "vcf_cache"   # per-gene het matrices cached here as .rds
OFFTARGET_FILE <- "results/offtarget/offtarget_summary.tsv"  # from 04_parse_offtargets.py


VCF_DIR      <- "data/1000G_highcov_GRCh38"
VCF_TEMPLATE <- file.path(
  VCF_DIR,
  "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)

dir.create(CACHE_DIR, showWarnings = FALSE)
dir.create("figures",  showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Load SNPs
# --------------------------------------------------------------------------- #

dat  <- read.delim(TSV_FILE, stringsAsFactors = FALSE, check.names = FALSE)
snp  <- dat[dat$variant_class == "common_snp" & dat$pass_basic == 1, ]
snp_u <- snp[!duplicated(snp$variant_id),
             c("gene", "variant_id", "het_rate", "ref", "alt")]
# variant_id is "1:pos:ref:alt" but VCFs use chr-prefixed contigs
snp_u$chrom_bare <- sub(":.*", "", snp_u$variant_id)
snp_u$chrom      <- paste0("chr", snp_u$chrom_bare)
snp_u$pos        <- as.integer(sapply(strsplit(snp_u$variant_id, ":"), `[[`, 2L))

cat(sprintf("Loaded %d unique common SNPs across %d genes\n",
            nrow(snp_u), length(unique(snp_u$gene))))

# Write sample IDs list for plot_ancestry_coverage.R (needs BCFtools in PATH)
SAMPLE_IDS_FILE <- file.path(CACHE_DIR, "sample_ids.txt")
if (!file.exists(SAMPLE_IDS_FILE)) {
  tmp_chrom <- snp_u$chrom[1]
  tmp_vcf   <- sprintf(VCF_TEMPLATE, sub("^chr", "", tmp_chrom))
  if (file.exists(tmp_vcf)) {
    ret <- system(sprintf("bcftools query -l '%s' > '%s'", tmp_vcf, SAMPLE_IDS_FILE))
    if (ret == 0) cat("Written:", SAMPLE_IDS_FILE, "\n")
  }
}

# --------------------------------------------------------------------------- #
# bcftools query: pull per-individual genotypes for a set of positions.
# Uses bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
# which emits one tab-delimited row per site: CHROM POS REF ALT GT1 GT2 ...
# Returns a logical matrix: rows = matched SNPs (pos:ref:alt), cols = individuals
# --------------------------------------------------------------------------- #

query_het_matrix <- function(chrom, pos_vec, ref_vec, alt_vec) {
  vcf <- sprintf(VCF_TEMPLATE, sub("^chr", "", chrom))  # template already has "chr"
  if (!file.exists(vcf)) stop("VCF not found: ", vcf)

  # Write a regions file for bcftools (chr:pos-pos, one per line)
  regions_file <- tempfile(fileext = ".regions")
  writeLines(paste0(chrom, "\t", pos_vec), regions_file)
  on.exit(unlink(regions_file), add = TRUE)

  cmd <- sprintf(
    "bcftools query -R '%s' -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' '%s' 2>/dev/null",
    regions_file, vcf
  )
  lines <- system(cmd, intern = TRUE)

  if (length(lines) == 0) {
    warning("No VCF lines returned for ", chrom); return(NULL)
  }

  rows <- lapply(lines, function(ln) {
    f     <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    pos_i <- as.integer(f[2])
    ref_i <- f[3]; alt_i <- f[4]
    gts   <- f[5:length(f)]
    het   <- gts %in% c("0|1","1|0","0/1","1/0")
    list(key = paste(pos_i, ref_i, alt_i, sep = ":"), het = het)
  })

  key_vcf <- sapply(rows, `[[`, "key")
  key_snp <- paste(pos_vec, ref_vec, alt_vec, sep = ":")
  matched <- match(key_snp, key_vcf)

  keep <- !is.na(matched)
  if (!any(keep)) { warning("No SNPs matched for ", chrom); return(NULL) }

  mat <- do.call(rbind, lapply(rows[matched[keep]], `[[`, "het"))
  rownames(mat) <- key_snp[keep]
  mat
}

# --------------------------------------------------------------------------- #
# Greedy set cover
# Returns data.frame with k (handles added) and cumulative coverage fraction
# --------------------------------------------------------------------------- #

greedy_coverage <- function(het_mat) {
  n_indiv  <- ncol(het_mat)
  covered  <- rep(FALSE, n_indiv)
  remaining <- seq_len(nrow(het_mat))
  result   <- data.frame(k = integer(0), frac_covered = numeric(0))

  while (length(remaining) > 0 && !all(covered)) {
    gains <- vapply(remaining,
                    function(i) sum(het_mat[i, ] & !covered),
                    integer(1))
    best     <- remaining[which.max(gains)]
    covered  <- covered | het_mat[best, ]
    result   <- rbind(result,
                      data.frame(k = nrow(result) + 1L,
                                 frac_covered = sum(covered) / n_indiv))
    remaining <- setdiff(remaining, best)
    # Stop early once gain is zero (all remaining SNPs redundant)
    if (max(gains) == 0) break
  }
  result
}

# --------------------------------------------------------------------------- #
# Main loop: per gene
# --------------------------------------------------------------------------- #

genes      <- sort(unique(snp_u$gene))
all_curves <- vector("list", length(genes))
names(all_curves) <- genes

for (g in genes) {
  cache_file <- file.path(CACHE_DIR, paste0(g, "_het_matrix.rds"))

  if (file.exists(cache_file)) {
    cat(sprintf("%-10s loading from cache\n", g))
    het_mat <- readRDS(cache_file)
  } else {
    cat(sprintf("%-10s querying 1KGP VCFs ...\n", g))
    gsnps  <- snp_u[snp_u$gene == g, ]
    chroms <- unique(gsnps$chrom)   # already chr-prefixed
    if (length(chroms) > 1) {
      cat("  WARNING: multiple chromosomes, taking first only\n")
      gsnps <- gsnps[gsnps$chrom == chroms[1], ]
    }
    het_mat <- query_het_matrix(gsnps$chrom[1],
                                gsnps$pos, gsnps$ref, gsnps$alt)
    if (is.null(het_mat)) {
      cat("  SKIPPING (no data)\n"); next
    }
    saveRDS(het_mat, cache_file)
    cat(sprintf("  %d SNPs matched, %d individuals\n",
                nrow(het_mat), ncol(het_mat)))
  }

  curve <- greedy_coverage(het_mat)
  curve$gene          <- g
  curve$n_snps_avail  <- nrow(het_mat)
  all_curves[[g]]     <- curve
}

coverage_df <- do.call(rbind, Filter(Negate(is.null), all_curves))
dir.create("gene_plot_inputs", showWarnings = FALSE)
write.csv(coverage_df, "figures/coverage_curves.csv",      row.names = FALSE)
write.csv(coverage_df, "gene_plot_inputs/coverage_curves.csv", row.names = FALSE)

# --------------------------------------------------------------------------- #
# Plot: one curve per gene
# --------------------------------------------------------------------------- #

# Mark the k at which each gene first exceeds 95%
marks_95 <- coverage_df %>%
  group_by(gene) %>%
  filter(frac_covered >= 0.95) %>%
  slice_min(k, n = 1) %>%
  ungroup()

p <- ggplot(coverage_df,
            aes(x = k, y = frac_covered * 100, colour = gene, group = gene)) +
  geom_hline(yintercept = 95, linetype = "dashed",
             colour = "grey55", linewidth = 0.4) +
  geom_line(linewidth = 0.75, alpha = 0.9) +
  geom_point(size = 1.4, alpha = 0.9) +
  geom_point(data = marks_95,
             aes(x = k, y = frac_covered * 100),
             shape = 21, fill = "white", size = 2.8, stroke = 0.9) +
  scale_y_continuous(limits = c(0, 101),
                     breaks = c(0, 25, 50, 75, 90, 95, 100),
                     labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = 1:max(coverage_df$k)) +
  labs(x     = "ASO targets added (greedily, to maximise coverage)",
       y     = "Individuals with ≥1 heterozygous ASO target",
       colour = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor  = element_blank(),
        legend.position   = "right",
        legend.text       = element_text(size = 8),
        legend.key.height = unit(0.5, "cm"))

cairo_pdf("figures/coverage_curves.pdf", width = 8, height = 5)
print(p)
dev.off()
cat("Written: figures/coverage_curves.pdf\n")

# Summary table: k needed for 80%, 90%, 95%, 99%
thresholds <- c(0.80, 0.90, 0.95, 0.99)
summ <- do.call(rbind, lapply(genes, function(g) {
  d <- coverage_df[coverage_df$gene == g, ]
  if (is.null(d) || nrow(d) == 0) return(NULL)
  ks <- sapply(thresholds, function(t) {
    w <- which(d$frac_covered >= t)
    if (length(w)) d$k[min(w)] else NA_integer_
  })
  setNames(data.frame(g, d$n_snps_avail[1], t(ks)),
           c("gene","n_snps", paste0("k_", thresholds*100, "pct")))
}))
cat("\n--- SNPs needed per gene ---\n")
print(summ, row.names = FALSE)

# --------------------------------------------------------------------------- #
# Coverage after specificity filtering
# Requires 04_parse_offtargets.py output with pass_specificity columns
# --------------------------------------------------------------------------- #

if (file.exists(OFFTARGET_FILE)) {
  ot <- read.delim(OFFTARGET_FILE, stringsAsFactors = FALSE)

  ot$pass_specificity_moderate <- as.logical(ot$pass_specificity_moderate)

  if (!"pass_specificity_moderate" %in% names(ot)) {
    cat("NOTE: offtarget_summary.tsv lacks pass_specificity columns; ",
        "re-run 04_parse_offtargets.py with --repeats to enable specificity-filtered coverage.\n")
  } else {
    # Moderate filter: no repeat + no exact (0mm) off-target + no paralog hits
    # A SNP passes if ANY allele passes across ALL search indices
    ot_snp <- ot[ot$variant_class == "common_snp", ]

    ot_per_allele <- aggregate(
      pass_specificity_moderate ~ variant_id + allele,
      data    = ot_snp,
      FUN     = all
    )
    ot_per_var <- aggregate(
      pass_specificity_moderate ~ variant_id,
      data    = ot_per_allele,
      FUN     = any
    )
    names(ot_per_var)[2] <- "pass_specificity"

    cat(sprintf("Specificity filter: %d/%d common SNPs pass moderate filter\n",
                sum(ot_per_var$pass_specificity), nrow(ot_per_var)))

    snp_u_spec <- merge(snp_u, ot_per_var, by = "variant_id", all.x = TRUE)
    snp_u_spec$pass_specificity[is.na(snp_u_spec$pass_specificity)] <- FALSE

    # Run greedy coverage for three SNP sets per gene
    run_coverage_set <- function(snp_subset, label) {
      out <- vector("list", length(genes))
      names(out) <- genes
      for (g in genes) {
        cache_file <- file.path(CACHE_DIR, paste0(g, "_het_matrix.rds"))
        if (!file.exists(cache_file)) next
        het_mat <- readRDS(cache_file)

        gsnps <- snp_subset[snp_subset$gene == g, ]
        if (nrow(gsnps) == 0) {
          out[[g]] <- data.frame(k = 0L, frac_covered = 0,
                                 gene = g, filter = label)
          next
        }
        # Subset het_mat rows to SNPs in gsnps that were in original matrix
        keep_rows <- rownames(het_mat) %in%
          paste(gsnps$pos, gsnps$ref, gsnps$alt, sep = ":")
        if (!any(keep_rows)) {
          out[[g]] <- data.frame(k = 0L, frac_covered = 0,
                                 gene = g, filter = label)
          next
        }
        curve      <- greedy_coverage(het_mat[keep_rows, , drop = FALSE])
        curve$gene <- g
        curve$filter <- label
        out[[g]] <- curve
      }
      do.call(rbind, Filter(Negate(is.null), out))
    }

    snp_raw  <- snp_u                              # no filter
    snp_basic <- snp_u                             # pass_basic already applied above
    snp_spec  <- snp_u_spec[snp_u_spec$pass_specificity == TRUE, ]

    cat("Running greedy coverage for three filter tiers...\n")
    df_raw   <- run_coverage_set(snp_raw,   "All common SNPs")
    df_basic <- run_coverage_set(snp_basic, "Pass basic filters")
    df_spec  <- run_coverage_set(snp_spec,  "Pass specificity filters")

    spec_all <- rbind(df_raw, df_basic, df_spec)
    spec_all$filter <- factor(spec_all$filter,
                              levels = c("All common SNPs",
                                         "Pass basic filters",
                                         "Pass specificity filters"))

    write.csv(spec_all, "figures/coverage_specificity.csv", row.names = FALSE)

    # ------------------------------------------------------------------ #
    # Per-superpopulation greedy coverage for specificity-filtered SNPs
    # ------------------------------------------------------------------ #
    PANEL_FILE   <- "integrated_call_samples_v3.20130502.ALL.panel"
    SAMPLE_ID_FILE <- file.path(CACHE_DIR, "sample_ids.txt")
    if (file.exists(PANEL_FILE) && file.exists(SAMPLE_ID_FILE)) {
      panel <- read.table(PANEL_FILE, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, fill = TRUE)
      vcf_samples <- readLines(SAMPLE_ID_FILE)   # ordered as in the VCF / het matrix cols

      superpops <- sort(unique(panel$super_pop))

      pop_spec_curves <- vector("list", length(genes) * length(superpops))
      idx_i <- 1L
      for (g in genes) {
        cache_file <- file.path(CACHE_DIR, paste0(g, "_het_matrix.rds"))
        if (!file.exists(cache_file)) next
        het_mat <- readRDS(cache_file)

        # Assign sample IDs to columns if missing
        if (is.null(colnames(het_mat)) && ncol(het_mat) == length(vcf_samples)) {
          colnames(het_mat) <- vcf_samples
        }

        gsnps <- snp_spec[snp_spec$gene == g, ]
        keep_rows <- rownames(het_mat) %in%
          paste(gsnps$pos, gsnps$ref, gsnps$alt, sep = ":")

        for (sp in superpops) {
          sp_samples <- panel$sample[panel$super_pop == sp]
          sp_cols    <- intersect(sp_samples, colnames(het_mat))
          if (length(sp_cols) == 0 || !any(keep_rows)) {
            pop_spec_curves[[idx_i]] <- data.frame(
              k = 0L, frac_covered = 0, gene = g, superpop = sp)
          } else {
            sub_mat <- het_mat[keep_rows, sp_cols, drop = FALSE]
            curve   <- greedy_coverage(sub_mat)
            curve$gene     <- g
            curve$superpop <- sp
            pop_spec_curves[[idx_i]] <- curve
          }
          idx_i <- idx_i + 1L
        }
      }

      pop_spec_df <- do.call(rbind, Filter(Negate(is.null), pop_spec_curves))
      pop_spec_df$gene     <- factor(pop_spec_df$gene,     levels = genes)
      pop_spec_df$superpop <- factor(pop_spec_df$superpop, levels = superpops)
      write.csv(pop_spec_df, "figures/pop_greedy_specificity.csv", row.names = FALSE)
      cat("Written: figures/pop_greedy_specificity.csv\n")
    }

    FILTER_COLS <- c(
      "All common SNPs"          = "#999999",
      "Pass basic filters"       = "#0072B2",
      "Pass specificity filters" = "#CC79A7"
    )

    pspec <- ggplot(spec_all,
                    aes(x = k, y = frac_covered * 100,
                        colour = filter, group = interaction(gene, filter))) +
      geom_hline(yintercept = 95, linetype = "dashed",
                 colour = "grey55", linewidth = 0.4) +
      geom_line(linewidth  = 0.65, alpha = 0.85) +
      geom_point(size = 1.2, alpha = 0.85) +
      scale_colour_manual(values = FILTER_COLS, name = NULL) +
      scale_y_continuous(limits = c(0, 101),
                         breaks = c(0, 25, 50, 75, 90, 95, 100),
                         labels = function(x) paste0(x, "%")) +
      scale_x_continuous(breaks = function(x) seq(1, max(x))) +
      labs(x      = "ASO targets added (greedy)",
           y      = "Population coverage",
           title  = "Coverage after specificity filtering") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            legend.position  = "bottom",
            legend.text      = element_text(size = 8))

    cairo_pdf("figures/coverage_specificity_filtered.pdf", width = 8, height = 5)
    print(pspec)
    dev.off()
    cat("Written: figures/coverage_specificity_filtered.pdf\n")
  }
} else {
  cat("NOTE: ", OFFTARGET_FILE, "not found; skipping specificity-filtered coverage.\n",
      "Run 03_search_offtargets.sl then 04_parse_offtargets.py first.\n")
}
