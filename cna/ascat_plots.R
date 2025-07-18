plot_ascat_allelic_segments <- function(segment_file, sample_id = NULL,
                                        exclude_chrXY = FALSE, min_seg_size = 0.5e6,
                                        offset = 0.07, line_width = 1.5,
                                        cn_cap = 5) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)

  # Define chromosome order
  chr_order <- c(as.character(1:22), "X", "Y")

  # Load and clean data
  segs <- read.delim(segment_file, stringsAsFactors = FALSE) %>%
    mutate(chr = toupper(as.character(chr)))

  if (exclude_chrXY) {
    segs <- segs %>% filter(!chr %in% c("X", "Y"))
  }

  if (!is.null(sample_id)) {
    segs <- segs %>% filter(sample == sample_id)
  }

  segs <- segs %>%
    mutate(seg_size = endpos - startpos) %>%
    filter(seg_size >= min_seg_size)

  segs <- segs %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    arrange(chr, startpos)

  # Create allele-specific long format with pair_id
  segs_with_pairs <- segs %>%
    mutate(pair_id = row_number()) %>%
    pivot_longer(cols = c(nMajor, nMinor), names_to = "allele", values_to = "copy_number") %>%
    mutate(chr = factor(chr, levels = chr_order))

  # Always offset alleles
  segs_with_pairs <- segs %>%
    mutate(pair_id = row_number()) %>%
    pivot_longer(cols = c(nMajor, nMinor), names_to = "allele", values_to = "copy_number") %>%
    mutate(chr = factor(chr, levels = chr_order),
           offset_cn = case_when(
             allele == "nMajor" ~ copy_number + offset,
             allele == "nMinor" ~ copy_number - offset,
             TRUE ~ copy_number
           )) %>%
    mutate(offset_cn = pmin(offset_cn, cn_cap))

  # Cap copy number at cn_cap
  segs_with_pairs <- segs_with_pairs %>%
    mutate(offset_cn = pmin(offset_cn, cn_cap))

  # Chromosome coordinates
  chr_lengths <- segs %>%
    group_by(chr) %>%
    summarise(chr_len = max(endpos), .groups = "drop") %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    arrange(chr) %>%
    mutate(chr_len = as.numeric(chr_len),
           chr_start = c(0, head(cumsum(chr_len), -1)))

  # Join coordinates
  segs_with_pairs <- segs_with_pairs %>%
    left_join(chr_lengths, by = "chr") %>%
    mutate(start_genome = startpos + chr_start,
           end_genome = endpos + chr_start)

  if (any(is.na(segs_with_pairs$chr_start))) {
    stop("NA values in 'chr_start'. Check chromosome naming.")
  }

  # Axis labels and vertical lines
  chr_midpoints <- chr_lengths %>%
    mutate(mid = chr_start + chr_len / 2)

  chr_boundaries <- chr_lengths %>%
    mutate(boundary = chr_start + chr_len) %>%
    pull(boundary)

  # Build custom y-axis with "5+" label
  y_breaks <- seq(0, cn_cap, by = 1)
  y_labels <- as.character(y_breaks)
  y_labels[length(y_labels)] <- paste0(cn_cap, "+")

  # Plot
  ggplot(segs_with_pairs) +
    geom_segment(aes(x = start_genome, xend = end_genome,
                     y = offset_cn, yend = offset_cn,
                     color = allele), linewidth = line_width) +
    geom_vline(xintercept = chr_boundaries, color = "gray90", linetype = "dashed", linewidth = 0.3) +
    scale_x_continuous(name = "Chromosome", breaks = chr_midpoints$mid, labels = chr_midpoints$chr) +
    scale_y_continuous(name = paste0("Allele-specific Copy Number"),
                       breaks = y_breaks, labels = y_labels, expand = c(0.01, 0), limits = c(-0.2, 5)) +
    scale_color_manual(values = c(nMajor = "#7D26CD", nMinor = "#00868B"), name = "Allele") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
}


setwd("/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/custom_ascat_plots")

pdf(file = "ascat_plot_TIL015_cfDNA.pdf", width = 7, height = 2.5)
plot_ascat_allelic_segments("/Volumes/datos_lab/ascat/TIL015/hg38_NeoPred/TIL015_cfDNA.segments.txt")
dev.off()

pdf(file = "ascat_plot_TIL015_FrTu.pdf", width = 7, height = 2.5)
plot_ascat_allelic_segments("/Volumes/datos_lab/ascat/TIL015/hg38_NeoPred/TIL015_FrTu.segments.txt")
dev.off()

pdf(file = "ascat_plot_Motricolor002_cfDNA.pdf", width = 7, height = 2.5)
plot_ascat_allelic_segments("/Volumes/datos_lab/ascat/Motricolor002/Motricolor002_cfDNA.segments.txt")
dev.off()

pdf(file = "ascat_plot_Motricolor002_FrTu.pdf", width = 7, height = 2.5)
plot_ascat_allelic_segments("/Volumes/datos_lab/ascat/Motricolor002/Motricolor002_FrTu.segments.txt")
dev.off()