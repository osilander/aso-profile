# Fig 1 — per-gene variant facets panel (STAT1/PIK3CA example locus)
# Output:
#   figures/adjusted_stat1_pik3ca_variant_facets_panel_swap.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(grid)
})

input_dir <- "gene_plot_inputs"
selected_genes <- c("STAT1", "PIK3CA")
output_pdf <- "figures/adjusted_stat1_pik3ca_variant_facets_panel_swap.pdf"

make_key <- function(gene, variant_class, pos, ref, alt) {
  paste(gene, variant_class, pos, ref, alt, sep = ":")
}

to_numeric <- function(x) suppressWarnings(as.numeric(x))

read_selected_contexts <- function(path) {
  x <- read.delim(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  names(x) <- paste0("V", seq_len(ncol(x)))
  
  x$VARIANT_CLASS <- x$V1
  x$GENE <- x$V2
  x$ID <- x$V3
  x$REF <- x$V4
  x$ALT <- x$V5
  x$MAF <- to_numeric(x$V12)
  x$HET_RATE <- to_numeric(x$V13)
  x$CONTEXT <- x$V14
  x$STRAND <- ifelse(x$VARIANT_CLASS == "common_snp",
                     x$V15,
                     sub(".*\\|([+-])::.*", "\\1", x$V29))
  x$GENOMIC_HGVS <- ifelse(x$VARIANT_CLASS == "common_snp", NA, x$V15)
  x$RNA_WINDOW <- x$V19
  x$POS <- ifelse(x$VARIANT_CLASS == "common_snp",
                  to_numeric(sub("^[^:]+:([0-9]+):.*", "\\1", x$ID)),
                  to_numeric(sub(".*:g\\.([0-9]+)[ACGT]>[ACGT].*", "\\1", x$GENOMIC_HGVS)))
  x$WINDOW_START <- to_numeric(sub(".*::chr[^:]+:([0-9]+)-[0-9]+\\([+-]\\).*", "\\1", x$V29))
  x$WINDOW_END <- to_numeric(sub(".*::chr[^:]+:[0-9]+-([0-9]+)\\([+-]\\).*", "\\1", x$V29))
  x$highlight_key <- make_key(x$GENE, x$VARIANT_CLASS, x$POS, x$REF, x$ALT)
  x
}

complement_base <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

in_intervals_by_gene <- function(gene, pos, features) {
  vapply(seq_along(pos), function(i) {
    f <- features[features$GENE == gene[i], ]
    nrow(f) > 0 && any(pos[i] > f$GENOMIC_START0 & pos[i] <= f$GENOMIC_END)
  }, logical(1))
}

cdna_label_for_row <- function(gene, pos, ref, alt, cds) {
  c <- cds[cds$GENE == gene, ]
  if (nrow(c) == 0) return(NA_character_)
  
  c <- c[order(c$TX_START), ]
  hit <- which(pos > c$GENOMIC_START0 & pos <= c$GENOMIC_END)
  if (length(hit) == 0) return(NA_character_)
  
  i <- hit[1]
  bases_before <- if (i == 1) 0 else sum(c$LENGTH[seq_len(i - 1)])
  
  if (c$STRAND[i] == "-") {
    offset <- c$GENOMIC_END[i] - pos
    ref <- complement_base(ref)
    alt <- complement_base(alt)
  } else {
    offset <- pos - (c$GENOMIC_START0[i] + 1)
  }
  
  paste0("c.", bases_before + offset + 1, ref, ">", alt)
}

assign_lanes <- function(pos, min_sep_bp) {
  lane <- rep(NA_integer_, length(pos))
  last_pos <- numeric(0)
  
  for (idx in order(pos)) {
    available <- which(pos[idx] - last_pos >= min_sep_bp)
    use_lane <- if (length(available) > 0) available[1] else length(last_pos) + 1
    
    if (use_lane > length(last_pos)) last_pos <- c(last_pos, -Inf)
    lane[idx] <- use_lane
    last_pos[use_lane] <- pos[idx]
  }
  
  lane
}

first_allele_sequence <- function(x) {
  no_brackets <- gsub("\\[([^]/]+)/[^]]+\\]", "\\1", x)
  gsub("[^ACGU]", "", no_brackets)
}

format_bp <- function(x) format(round(x), big.mark = ",", scientific = FALSE)

wrap_label <- function(x, width = 18) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

dat <- read.delim(file.path(input_dir, "all_genes.variants_for_plot.tsv"),
                  stringsAsFactors = FALSE,
                  check.names = FALSE)
exons <- read.delim(file.path(input_dir, "all_genes.exons_for_plot.tsv"),
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
cds <- read.delim(file.path(input_dir, "all_genes.cds_for_plot.tsv"),
                  stringsAsFactors = FALSE,
                  check.names = FALSE)
tx <- read.delim(file.path(input_dir, "all_genes.transcripts.tsv"),
                 stringsAsFactors = FALSE,
                 check.names = FALSE)
summary_counts <- read.delim(file.path(input_dir, "all_genes.summary.tsv"),
                             stringsAsFactors = FALSE,
                             check.names = FALSE)
selected_contexts <- read_selected_contexts("selected_STAT1_PIK3CA_contexts.tsv")

for (col in intersect(c("POS", "START0", "END", "MAF", "HET_RATE", "TRANSCRIPT_SPAN_BP"), names(dat))) {
  dat[[col]] <- to_numeric(dat[[col]])
}
for (col in intersect(c("GENOMIC_START0", "GENOMIC_END", "TX_START", "TX_END"), names(exons))) {
  exons[[col]] <- to_numeric(exons[[col]])
}
for (col in intersect(c("GENOMIC_START0", "GENOMIC_END", "TX_START", "TX_END", "LENGTH"), names(cds))) {
  cds[[col]] <- to_numeric(cds[[col]])
}
for (col in c("TRANSCRIPT_START0", "TRANSCRIPT_END", "TRANSCRIPT_SPAN_BP")) {
  tx[[col]] <- to_numeric(tx[[col]])
}

genes <- selected_genes
dat <- dat[dat$GENE %in% genes, ]
exons <- exons[exons$GENE %in% genes, ]
cds <- cds[cds$GENE %in% genes, ]
tx <- tx[tx$GENE %in% genes, ]
selected_contexts <- selected_contexts[selected_contexts$GENE %in% genes, ]

dat$highlight_key <- make_key(dat$GENE, dat$VARIANT_CLASS, dat$POS, dat$REF, dat$ALT)

# Add selected context variants that are outside the selected transcript table.
missing_selected <- selected_contexts[!selected_contexts$highlight_key %in% dat$highlight_key, ]
if (nrow(missing_selected) > 0) {
  extra <- as.data.frame(setNames(lapply(names(dat), function(x) rep(NA, nrow(missing_selected))),
                                  names(dat)),
                         stringsAsFactors = FALSE)
  extra$GENE <- missing_selected$GENE
  extra$VARIANT_CLASS <- missing_selected$VARIANT_CLASS
  extra$CHROM <- sub("^(chr[^:]+):.*", "\\1",
                     sub(".*::(chr[^:]+:[0-9]+-[0-9]+\\([+-]\\)).*", "\\1", missing_selected$V29))
  extra$POS <- missing_selected$POS
  extra$START0 <- missing_selected$POS - 1
  extra$END <- missing_selected$POS
  extra$REF <- missing_selected$REF
  extra$ALT <- missing_selected$ALT
  extra$ID <- missing_selected$ID
  extra$RS_ID <- NA
  extra$LABEL <- missing_selected$ID
  extra$CONTEXT <- missing_selected$CONTEXT
  extra$MAF <- missing_selected$MAF
  extra$HET_RATE <- missing_selected$HET_RATE
  extra$GENOMIC_HGVS <- missing_selected$GENOMIC_HGVS
  extra$TRANSCRIPT_ID <- tx$TRANSCRIPT_ID[match(missing_selected$GENE, tx$GENE)]
  extra$TRANSCRIPT_NAME <- tx$TRANSCRIPT_NAME[match(missing_selected$GENE, tx$GENE)]
  extra$TRANSCRIPT_SPAN_BP <- tx$TRANSCRIPT_SPAN_BP[match(missing_selected$GENE, tx$GENE)]
  extra$STRAND <- missing_selected$STRAND
  extra$highlight_key <- missing_selected$highlight_key
  dat <- rbind(dat, extra)
}

for (obj in c("dat", "exons", "cds", "tx")) {
  tmp <- get(obj)
  tmp$GENE <- factor(tmp$GENE, levels = genes)
  assign(obj, tmp)
}

dat$plot_context <- ifelse(in_intervals_by_gene(dat$GENE, dat$POS, cds),
                           "CDS",
                           ifelse(in_intervals_by_gene(dat$GENE, dat$POS, exons),
                                  "exon_nonCDS",
                                  "intron"))
dat$plot_group <- ifelse(dat$VARIANT_CLASS == "pathogenic_missense",
                         "P/LP missense",
                         "Common SNP")
dat$location_shape <- ifelse(dat$plot_context == "intron", "Intronic", "Exonic")
dat$variant_label <- paste0(format_bp(dat$POS), " ", dat$REF, ">", dat$ALT)
dat$cdna_label <- vapply(seq_len(nrow(dat)), function(i) {
  cdna_label_for_row(dat$GENE[i], dat$POS[i], dat$REF[i], dat$ALT[i], cds)
}, character(1))
dat$path_label <- ifelse(!is.na(dat$cdna_label), dat$cdna_label, dat$variant_label)
dat$plot_label <- ifelse(dat$plot_group == "Common SNP" &
                           !is.na(dat$RS_ID) &
                           dat$RS_ID != "NA" &
                           dat$RS_ID != "",
                         dat$RS_ID,
                         dat$variant_label)
dat$highlight_key <- make_key(dat$GENE, dat$VARIANT_CLASS, dat$POS, dat$REF, dat$ALT)
dat$highlight_target <- dat$highlight_key %in% selected_contexts$highlight_key

lane_dat <- do.call(rbind, lapply(split(dat, list(dat$GENE, dat$plot_group), drop = TRUE), function(x) {
  span <- max(x$POS, na.rm = TRUE) - min(x$POS, na.rm = TRUE)
  min_sep <- max(800, span * 0.018)
  x$lane <- assign_lanes(x$POS, min_sep)
  x$y <- ifelse(x$plot_group == "Common SNP",
                0.18 + (x$lane - 1) * 0.055,
                -0.18 - (x$lane - 1) * 0.055)
  x
}))

dat <- lane_dat[order(lane_dat$GENE, lane_dat$POS), ]
common <- dat[dat$plot_group == "Common SNP", ]
path <- dat[dat$plot_group == "P/LP missense", ]

common$label_me <- common$highlight_target
path$label_me <- FALSE
for (g in levels(dat$GENE)) {
  idx <- which(path$GENE == g)
  if (length(idx) == 0) next
  
  idx <- idx[order(path$POS[idx])]
  cluster_id <- cumsum(c(TRUE, diff(path$POS[idx]) > 100))
  for (cl in unique(cluster_id)) {
    cluster_idx <- idx[cluster_id == cl]
    if (length(cluster_idx) > 2) {
      keep_n <- max(1, ceiling(length(cluster_idx) / 4))
      keep <- cluster_idx[round(seq(1, length(cluster_idx), length.out = keep_n))]
    } else {
      keep <- cluster_idx
    }
    path$label_me[keep] <- TRUE
  }
}
path$label_me[path$highlight_target] <- TRUE

label_common <- common[common$label_me, ]
highlight_variants <- dat[dat$highlight_target, ]

exons$xmin <- exons$GENOMIC_START0 + 1
exons$xmax <- exons$GENOMIC_END
tx$xmin <- tx$TRANSCRIPT_START0 + 1
tx$xmax <- tx$TRANSCRIPT_END

sc <- summary_counts[summary_counts$GENE %in% genes, ]
facet_labels <- setNames(
  paste0(sc$GENE, "  pathogenic/likely pathogenic SNVs: ", sc$PATHOGENIC_IN_SELECTED_TRANSCRIPT,
         "; common SNPs: ", sc$COMMON_IN_SELECTED_TRANSCRIPT),
  sc$GENE
)

selected_contexts$gene_order <- match(selected_contexts$GENE, genes)
selected_contexts$panel_order <- ifelse(selected_contexts$VARIANT_CLASS == "pathogenic_missense", 1, 2)
selected_contexts <- selected_contexts[order(selected_contexts$gene_order,
                                             selected_contexts$panel_order,
                                             selected_contexts$POS), ]
selected_contexts$Type <- ifelse(selected_contexts$VARIANT_CLASS == "common_snp", "SNP", "P/LP")
selected_contexts$plot_path_label <- dat$path_label[match(selected_contexts$highlight_key, dat$highlight_key)]
selected_contexts$ASO <- paste0("ASO", ave(selected_contexts$POS, selected_contexts$GENE, FUN = seq_along))
selected_contexts$Variant <- ifelse(
  selected_contexts$VARIANT_CLASS == "common_snp",
  paste0("SNP ", format_bp(selected_contexts$POS), " ", selected_contexts$REF, ">", selected_contexts$ALT),
  paste0("P/LP ", selected_contexts$plot_path_label)
)
selected_contexts$MAF_label <- ifelse(is.na(selected_contexts$MAF), "NA", sprintf("%.3f", selected_contexts$MAF))
selected_contexts$RNA_WINDOW_BOLD <- sub("(\\[[^]]+\\])", "<b>\\1</b>", selected_contexts$RNA_WINDOW)
selected_contexts$RNA_FIRST_ALLELE <- first_allele_sequence(selected_contexts$RNA_WINDOW)
selected_contexts$GC_percent <- 100 * nchar(gsub("[^GC]", "", selected_contexts$RNA_FIRST_ALLELE)) /
  nchar(selected_contexts$RNA_FIRST_ALLELE)
selected_contexts$GC_percent_label <- sprintf("%.0f%%", selected_contexts$GC_percent)
selected_contexts$GC_status <- ifelse(selected_contexts$GC_percent >= 40 & selected_contexts$GC_percent <= 60,
                                      "OK",
                                      "X")
selected_contexts$GC_col <- ifelse(selected_contexts$GC_status == "OK", "#238B45", "#CB181D")

table_panel <- data.frame(
  ASO = selected_contexts$ASO,
  Gene = selected_contexts$GENE,
  Variant = selected_contexts$Variant,
  Type = selected_contexts$Type,
  MAF = selected_contexts$MAF_label,
  Context = selected_contexts$RNA_WINDOW_BOLD,
  GC = selected_contexts$GC_percent_label,
  GC_status = selected_contexts$GC_status,
  GC_col = selected_contexts$GC_col,
  y = rev(seq_len(nrow(selected_contexts))),
  stringsAsFactors = FALSE
)
table_header <- data.frame(
  label = c("ASO", "Gene", "Variant", "MAF", "RNA 17-mer context", "GC"),
  x = c(0.02, 0.13, 0.27, 0.55, 0.61, 0.92),
  hjust = c(0, 0, 0, 0, 0, 1),
  stringsAsFactors = FALSE
)
table_cells <- rbind(
  data.frame(x = 0.02, y = table_panel$y, label = table_panel$ASO, hjust = 0),
  data.frame(x = 0.13, y = table_panel$y, label = table_panel$Gene, hjust = 0),
  data.frame(x = 0.27, y = table_panel$y, label = table_panel$Variant, hjust = 0),
  data.frame(x = 0.55, y = table_panel$y, label = table_panel$MAF, hjust = 0),
  data.frame(x = 0.92, y = table_panel$y, label = table_panel$GC, hjust = 1)
)
table_context_cells <- data.frame(x = 0.61, y = table_panel$y, label = table_panel$Context, hjust = 0)
table_gc_cells <- data.frame(x = 0.955,
                             y = table_panel$y,
                             label = table_panel$GC_status,
                             colour = table_panel$GC_col,
                             hjust = 0)
hap_variants <- selected_contexts
hap_variants$GENE <- factor(hap_variants$GENE, levels = genes)
hap_variants$Type <- ifelse(hap_variants$VARIANT_CLASS == "common_snp", "SNP", "P/LP")
hap_variants$point_y <- 0
hap_variants$aso_id <- ave(hap_variants$POS, hap_variants$GENE, FUN = seq_along)
hap_variants$aso_lane <- unlist(lapply(split(hap_variants$POS, hap_variants$GENE), assign_lanes, min_sep_bp = 1200),
                                use.names = FALSE)
hap_variants$aso_y <- 0.070 + (hap_variants$aso_lane - 1) * 0.060
hap_variants$aso_label <- paste0("ASO", hap_variants$aso_id)
hap_variants$ASO_START <- hap_variants$POS - 500
hap_variants$ASO_END <- hap_variants$POS + 500
hap_variants$aso_colour <- ifelse(hap_variants$Type == "SNP", "#0072B2", "grey35")
hap_variants$aso_label_y <- hap_variants$aso_y + 0.075
hap_variants$aso_label_y[hap_variants$GENE == "PIK3CA" & hap_variants$aso_label == "ASO3"] <- 0.215
hap_variants$aso_label_y[hap_variants$GENE == "PIK3CA" & hap_variants$aso_label == "ASO4"] <- 0.310

hap_exons <- exons
hap_exons$ymin <- -0.04
hap_exons$ymax <- 0.04

hap_tx <- tx
hap_tx$y <- 0

hap_facet_labels <- setNames(paste0(genes, " patient haplotype"), genes)

y_min <- min(-0.90, min(dat$y, na.rm = TRUE) - 0.25)
y_max <- max(0.90, max(dat$y, na.rm = TRUE) + 0.25)

group_text <- do.call(rbind, lapply(split(tx, tx$GENE), function(x) {
  span <- x$xmax[1] - x$xmin[1]
  data.frame(GENE = x$GENE[1],
             x = x$xmin[1] + 0.025 * span,
             y = c(0.72, -0.48),
             label = c("common SNPs", "pathogenic SNVs"),
             colour = c("#E6550D", "grey45"),
             stringsAsFactors = FALSE)
}))
group_text$GENE <- factor(group_text$GENE, levels = genes)

legend_gene <- "PIK3CA"
legend_tx <- tx[tx$GENE == legend_gene, ][1, ]
legend_x0 <- legend_tx$xmin + 0.025 * (legend_tx$xmax - legend_tx$xmin)
legend_x1 <- legend_tx$xmin + 0.18 * (legend_tx$xmax - legend_tx$xmin)
legend_y0 <- y_min + 0.19
legend_y1 <- legend_y0 + 0.13
maf_cols <- viridisLite::inferno(60, begin = 0.12, end = 0.95, direction = -1)
legend_breaks <- seq(legend_x0, legend_x1, length.out = length(maf_cols) + 1)
maf_legend <- data.frame(
  GENE = factor(rep(legend_gene, length(maf_cols)), levels = genes),
  xmin = legend_breaks[-length(legend_breaks)],
  xmax = legend_breaks[-1],
  ymin = legend_y0,
  ymax = legend_y1,
  fill = maf_cols
)
maf_legend_text <- data.frame(
  GENE = factor(rep(legend_gene, 3), levels = genes),
  x = c(legend_x0, legend_x1, (legend_x0 + legend_x1) / 2),
  y = c(legend_y0 - 0.075, legend_y0 - 0.075, legend_y1 + 0.075),
  label = c("0.2", "0.5", "MAF"),
  hjust = c(0, 1, 0.5)
)

flow_steps <- data.frame(
  x = seq(0.075, 0.925, length.out = 5),
  y = 0.50,
  label = wrap_label(c("Patient WGS",
                       "Phase pathogenic variant",
                       "Identify phased heterozygous transcribed SNPs",
                       "Match prevalidated ASO to patient phased SNPs",
                       "Test allele-specific knockdown in patient cells"),
                     width = 16),
  stringsAsFactors = FALSE
)
flow_box_width <- 0.17
flow_box_height <- 0.62
flow_steps$xmin <- flow_steps$x - flow_box_width / 2
flow_steps$xmax <- flow_steps$x + flow_box_width / 2
flow_steps$ymin <- flow_steps$y - flow_box_height / 2
flow_steps$ymax <- flow_steps$y + flow_box_height / 2
flow_arrows <- data.frame(
  x = head(flow_steps$xmax, -1) + 0.004,
  xend = tail(flow_steps$xmin, -1) - 0.004,
  y = 0.50,
  yend = 0.50
)

draw_flowchart <- function() {
  grid.segments(x0 = unit(flow_arrows$x, "npc"),
                x1 = unit(flow_arrows$xend, "npc"),
                y0 = unit(flow_arrows$y, "npc"),
                y1 = unit(flow_arrows$yend, "npc"),
                arrow = arrow(length = unit(7, "pt"), type = "closed"),
                gp = gpar(col = "grey35", lwd = 0.80, lineend = "round"))
  for (i in seq_len(nrow(flow_steps))) {
    grid.roundrect(x = unit(flow_steps$x[i], "npc"),
                   y = unit(flow_steps$y[i], "npc"),
                   width = unit(flow_box_width, "npc"),
                   height = unit(flow_box_height, "npc"),
                   r = unit(5, "pt"),
                   gp = gpar(fill = "grey96", col = "grey35", lwd = 0.22))
    grid.text(flow_steps$label[i],
              x = unit(flow_steps$x[i], "npc"),
              y = unit(flow_steps$y[i], "npc"),
              gp = gpar(col = "grey18", fontsize = 11.5, lineheight = 0.88))
  }
}

p <- ggplot() +
  geom_segment(data = tx,
               aes(x = xmin, xend = xmax, y = 0, yend = 0),
               colour = "grey35",
               linewidth = 0.35) +
  geom_rect(data = exons,
            aes(xmin = xmin, xmax = xmax, ymin = -0.045, ymax = 0.045),
            fill = "grey70",
            colour = "grey30",
            linewidth = 0.18) +
  geom_point(data = common,
             aes(x = POS, y = y, fill = MAF, shape = location_shape),
             colour = "grey18",
             size = 2.3,
             stroke = 0.32) +
  geom_point(data = path,
             aes(x = POS, y = y, shape = location_shape),
             fill = "grey60",
             colour = "grey25",
             size = 2.1,
             stroke = 0.29,
             alpha = 0.75) +
  geom_point(data = highlight_variants,
             aes(x = POS, y = y),
             shape = 21,
             fill = NA,
             colour = "#00A6D6",
             size = 3.5,
             stroke = 1.25,
             show.legend = FALSE) +
  geom_rect(data = maf_legend,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = maf_legend$fill,
            colour = NA,
            inherit.aes = FALSE) +
  geom_rect(data = data.frame(GENE = factor(legend_gene, levels = genes),
                              xmin = legend_x0,
                              xmax = legend_x1,
                              ymin = legend_y0,
                              ymax = legend_y1),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA,
            colour = "grey25",
            linewidth = 0.15,
            inherit.aes = FALSE) +
  geom_text(data = maf_legend_text,
            aes(x = x, y = y, label = label, hjust = hjust),
            size = 3.05,
            colour = "grey20",
            inherit.aes = FALSE) +
  geom_text(data = group_text,
            aes(x = x, y = y, label = label, colour = colour),
            hjust = 0,
            fontface = "bold",
            size = 4.3,
            inherit.aes = FALSE) +
  geom_text_repel(data = label_common,
                  aes(x = POS, y = y, label = plot_label),
                  colour = "#7F2704",
                  size = 2.65,
                  nudge_y = 0.22,
                  min.segment.length = 0,
                  segment.color = "#7F2704",
                  segment.size = 0.12,
                  segment.alpha = 0.75,
                  box.padding = 0.40,
                  point.padding = 0.35,
                  force = 10,
                  force_pull = 0.01,
                  max.time = 8,
                  max.iter = 30000,
                  max.overlaps = Inf,
                  seed = 3) +
  facet_wrap(~ GENE, scales = "free_x", ncol = 1, labeller = labeller(GENE = facet_labels)) +
  scale_x_continuous(labels = format_bp, expand = expansion(mult = c(0.015, 0.015))) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = NULL, expand = expansion(mult = 0)) +
  scale_fill_viridis_c(option = "inferno",
                       begin = 0.12,
                       end = 0.95,
                       direction = -1,
                       na.value = "#D95F02",
                       guide = "none") +
  scale_shape_manual(values = c("Intronic" = 21, "Exonic" = 24), guide = "none") +
  scale_colour_identity() +
  labs(x = "Genomic position (bp)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = "grey70", linewidth = 0.25),
        strip.text = element_text(face = "plain", size = 12.6),
        legend.position = "none",
        panel.spacing = unit(0.9, "lines"),
        plot.margin = margin(8, 8, 12, 42))

hap_plot <- ggplot() +
  geom_segment(data = hap_tx,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               colour = "grey35",
               linewidth = 0.35) +
  geom_rect(data = hap_exons,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey70",
            colour = "grey30",
            linewidth = 0.18) +
  geom_segment(data = hap_variants,
               aes(x = ASO_START, xend = ASO_END, y = aso_y, yend = aso_y, colour = aso_colour),
               linewidth = 1.15,
               lineend = "round") +
  geom_segment(data = hap_variants,
               aes(x = POS, xend = POS, y = aso_y, yend = aso_label_y - 0.025),
               colour = "grey45",
               linewidth = 0.08,
               alpha = 0.70) +
  geom_point(data = hap_variants[hap_variants$Type == "SNP", ],
             aes(x = POS, y = point_y, fill = MAF),
             shape = 21,
             colour = "grey18",
             size = 3.35,
             stroke = 0.35) +
  geom_point(data = hap_variants[hap_variants$Type == "P/LP", ],
             aes(x = POS, y = point_y),
             shape = 24,
             fill = "grey60",
             colour = "grey25",
             size = 2.95,
             stroke = 0.46) +
  geom_text(data = hap_variants,
            aes(x = POS, y = aso_label_y, label = aso_label),
            size = 2.95,
            colour = "grey20") +
  facet_wrap(~ GENE, scales = "free_x", ncol = 1, labeller = labeller(GENE = hap_facet_labels)) +
  scale_x_continuous(labels = format_bp, expand = expansion(mult = c(0.015, 0.015))) +
  scale_y_continuous(limits = c(-0.10, 0.38), breaks = NULL, expand = expansion(mult = 0)) +
  scale_fill_viridis_c(option = "inferno",
                       begin = 0.12,
                       end = 0.95,
                       direction = -1,
                       limits = c(0.2, 0.5),
                       na.value = "#D95F02",
                       guide = "none") +
  scale_colour_identity() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = "grey75", linewidth = 0.2),
        strip.text = element_text(face = "plain", size = 11.3),
        panel.spacing = unit(0.70, "lines"),
        plot.margin = margin(10, 8, 12, 42))

table_plot <- ggplot() +
  geom_rect(data = table_panel,
            aes(xmin = 0, xmax = 1, ymin = y - 0.42, ymax = y + 0.42, fill = Type),
            alpha = 0.08,
            colour = NA) +
  geom_hline(yintercept = table_panel$y + 0.5, colour = "grey88", linewidth = 0.2) +
  geom_text(data = table_header,
            aes(x = x, y = nrow(table_panel) + 1, label = label, hjust = hjust),
            fontface = "bold",
            size = 3.95,
            colour = "grey15") +
  geom_text(data = table_cells,
            aes(x = x, y = y, label = label, hjust = hjust),
            size = 3.65,
            colour = "grey20") +
  geom_richtext(data = table_context_cells,
                aes(x = x, y = y, label = label, hjust = hjust),
                size = 3.65,
                colour = "grey20",
                fill = NA,
                label.color = NA,
                label.padding = unit(0, "pt")) +
  geom_text(data = table_gc_cells,
            aes(x = x, y = y, label = label, hjust = hjust, colour = colour),
            size = 3.95,
            fontface = "bold") +
  scale_colour_identity() +
  scale_fill_manual(values = c("P/LP" = "grey40", "SNP" = "#D95F02"), guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.3, nrow(table_panel) + 1.5), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(12, 16, 8, 42))

cairo_pdf(output_pdf, width = 10, height = 13.8)
grid.newpage()
print(p, vp = viewport(x = 0.5, y = 1, width = 1, height = 0.425, just = c("center", "top")))
print(hap_plot, vp = viewport(x = 0.5, y = 0.365, width = 1, height = 0.200, just = c("center", "bottom")))
print(table_plot, vp = viewport(x = 0.5, y = 0.135, width = 1, height = 0.215, just = c("center", "bottom")))
pushViewport(viewport(x = 0.5, y = 0, width = 1, height = 0.125, just = c("center", "bottom")))
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.90, height = 1, just = c("center", "center")))
draw_flowchart()
popViewport(2)
grid.text("a", x = unit(0.018, "npc"), y = unit(0.968, "npc"),
          just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 21, col = "grey15"))
grid.text("b", x = unit(0.018, "npc"), y = unit(0.770, "npc"),
          just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 21, col = "grey15"))
grid.text("c", x = unit(0.018, "npc"), y = unit(0.565, "npc"),
          just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 21, col = "grey15"))
grid.text("d", x = unit(0.018, "npc"), y = unit(0.350, "npc"),
          just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 21, col = "grey15"))
grid.text("e", x = unit(0.018, "npc"), y = unit(0.125, "npc"),
          just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 21, col = "grey15"))
dev.off()
