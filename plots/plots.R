# CCF compare plots:
# compare CCF data from two different tumor samples (e.g. biopsy vs cfDNA)

ccf_compare_plot <- function(
    ccf_tsv,
    patient,
    dna_presence_criteria = "vaf", # vaf or alt_reads
    min_vaf = 0.01, # minimum VAF to consider a mutation present in a sample
    min_alt_reads = 1, # minimum ALT reads to consider a mutation present in a sample
    show_immunogenic = TRUE, # whether to highlight immunogenic mutations
    label_immunogenic = TRUE, # whether to label immunogenic mutations
    limit_ccf = TRUE # whether to limit CCF values to max 1
    ) {
    # load libraries
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(ggridges))
    suppressPackageStartupMessages(library(ggrepel))
    suppressPackageStartupMessages(library(patchwork))

    # load ccf tsv
    ccf <- fread(ccf_tsv)

    # Define thresholds depending on criteria
    if (dna_presence_criteria == "vaf") {
        present_FrTu <- ccf$vaf_Tumor_1 >= min_vaf
        present_cfDNA <- ccf$vaf_Tumor_2 >= min_vaf
    } else if (dna_presence_criteria == "alt_reads") {
        present_FrTu <- ccf$ALT_counts_FrTu >= min_alt_reads
        present_cfDNA <- ccf$ALT_counts_cfDNA >= min_alt_reads
    } else {
        stop("dna_presence_criteria must be either 'vaf' or 'alt_reads'")
    }

    # Add dna_presence classification
    ccf <- ccf %>%
        mutate(dna_presence = case_when(
            present_FrTu & present_cfDNA ~ "Shared",
            present_FrTu & !present_cfDNA ~ "FrTu_only",
            !present_FrTu & present_cfDNA ~ "cfDNA_only",
            TRUE ~ "undetermined"
        ))

    # Handle limit_ccf
    if (limit_ccf) {
        ccf <- ccf %>%
            mutate(
                ccf_FrTu  = pmin(ccf_FrTu, 1),
                ccf_cfDNA = pmin(ccf_cfDNA, 1)
            )
    }

    # Base scatter plot
    p <- ggplot(ccf, aes(x = ccf_FrTu, y = ccf_cfDNA, color = dna_presence)) +
        geom_point(alpha = 0.3) +
        scale_color_manual(values = c("cfDNA_only" = "#FCC557", "FrTu_only" = "#D9D8D4", "Shared" = "#7CC6BF")) +
        theme_classic() +
        theme(legend.position = "none")

    # Conditionally add immunogenic labels
    if (show_immunogenic && "immunogenic" %in% colnames(ccf)) {
        p <- p +
            geom_point(
                data = ccf %>% filter(immunogenic == "yes"),
                aes(x = ccf_FrTu, y = ccf_cfDNA),
                shape = 8, color = "#333333", size = 2, inherit.aes = FALSE
            )
    }

    if (show_immunogenic && label_immunogenic && "immunogenic" %in% colnames(ccf)) {
        p <- p +
            geom_text_repel(
                data = ccf %>% filter(immunogenic == "yes"),
                aes(x = ccf_FrTu, y = ccf_cfDNA, label = gene_name),
                color = "grey30", size = 3.3,
                min.segment.length = unit(0, "lines"),
                segment.color = "grey",
                max.overlaps = Inf,
                box.padding = 0.5,
                fontface = "italic",
                force = 100,
                max.time = 50,
                inherit.aes = FALSE
            )
    }

    # Apply axis labels and limits
    if (limit_ccf) {
        p <- p +
            scale_x_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25)) +
            scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25)) +
            xlab("CCF FrTu") +
            ylab("CCF cfDNA")
    } else {
        p <- p +
            xlab("Variant Copy Number FrTu") +
            ylab("Variant Copy Number cfDNA")
    }

    # Density and histogram plots also respect limits
    dens1 <- ggplot(ccf, aes(x = ccf_FrTu, y = dna_presence, fill = dna_presence)) +
        geom_density_ridges(alpha = 0.4) +
        scale_fill_manual(values = c("#FCC557", "#D9D8D4", "#7CC6BF")) +
        theme_void() +
        theme(legend.position = "none")

    if (limit_ccf) {
        dens1 <- dens1 + xlim(-0.05, 1.05)
    }

    dens2 <- ggplot(ccf, aes(x = ccf_cfDNA, y = dna_presence, fill = dna_presence)) +
        geom_density_ridges(alpha = 0.4) +
        scale_fill_manual(values = c("#FCC557", "#D9D8D4", "#7CC6BF")) +
        theme_void() +
        # theme(legend.position = "none") +
        coord_flip()

    if (limit_ccf) {
        dens2 <- dens2 + xlim(-0.05, 1.05)
    }

    hist1 <- ggplot(ccf, aes(x = ccf_FrTu, fill = dna_presence)) +
        geom_histogram(binwidth = 0.05) +
        scale_fill_manual(values = c("#FCC557", "#D9D8D4", "#7CC6BF")) +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            legend.position = "none"
        ) +
        scale_y_reverse()

    if (limit_ccf) {
        hist1 <- hist1 + xlim(-0.05, 1.05)
    }

    hist2 <- ggplot(ccf, aes(x = ccf_cfDNA, fill = dna_presence)) +
        geom_histogram(binwidth = 0.05) +
        scale_fill_manual(values = c("#FCC557", "#D9D8D4", "#7CC6BF")) +
        theme_classic() +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.x = element_text(angle = 90),
            legend.position = "none"
        ) +
        coord_flip() +
        scale_y_reverse()

    if (limit_ccf) {
        hist2 <- hist2 + xlim(-0.05, 1.05)
    }

    # Combine panels
    patch_plot <- plot_spacer() + dens1 + plot_spacer() +
        hist2 + p + dens2 +
        plot_spacer() + hist1 + plot_spacer() +
        plot_layout(ncol = 3, nrow = 3, widths = c(1, 6.5, 1), heights = c(1, 6.5, 1)) +
        plot_annotation(title = patient)

    print(patch_plot)
}


# ASCAT segment plots:
# customizable function to plot ASCAT allele-specific segments

plot_allelic_segments <- function(
    segment_file,
    nmaj_color = "#7D26CD",
    nmin_color = "#00868B",
    sample_id = NULL,
    exclude_chrXY = FALSE,
    min_seg_size = 0.5e6,
    offset = 0.07,
    line_width = 1.5,
    cn_cap = 5) {
    # load libraries
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(scales))
    suppressPackageStartupMessages(library(data.table))

    # Define chromosome order
    chr_order <- c(as.character(1:22), "X", "Y")

    # Load data (robust to large TSVs)
    segs <- suppressWarnings(as_tibble(fread(segment_file)))

    # Detect format (ASCAT vs PURPLE) and normalize columns:
    # Expected unified columns after normalization: chr, startpos, endpos, nMajor, nMinor, sample
    nm <- names(segs)
    if ((("chromosome" %in% nm) || ("Chromosome" %in% nm)) &&
        ("majorAlleleCopyNumber" %in% nm) && ("minorAlleleCopyNumber" %in% nm)) {
        # PURPLE format
        chrom_col <- if ("chromosome" %in% nm) "chromosome" else "Chromosome"
        # prefer canonical names, with fallbacks sometimes seen in outputs
        start_col <- if ("start" %in% nm) {
            "start"
        } else if ("Start" %in% nm) {
            "Start"
        } else if ("minStart" %in% nm) {
            "minStart"
        } else {
            stop("Could not find a start column in PURPLE file (expected one of: start, Start, minStart)")
        }
        end_col <- if ("end" %in% nm) {
            "end"
        } else if ("End" %in% nm) {
            "End"
        } else if ("maxStart" %in% nm) {
            # rare case in SV-annotated rows; use maxStart as end
            "maxStart"
        } else {
            stop("Could not find an end column in PURPLE file (expected one of: end, End, maxStart)")
        }

        segs <- segs %>%
            mutate(
                chr = toupper(gsub("^CHR", "", toupper(as.character(.data[[chrom_col]])))),
                startpos = as.numeric(.data[[start_col]]),
                endpos = as.numeric(.data[[end_col]]),
                nMajor = as.numeric(round(majorAlleleCopyNumber)),
                nMinor = as.numeric(round(minorAlleleCopyNumber))
            )
    } else {
        # ASCAT-like: try to map common variants
        if (!("chr" %in% nm) && ("chromosome" %in% nm)) segs <- segs %>% rename(chr = chromosome)
        if (!("chr" %in% names(segs)) && ("Chromosome" %in% nm)) segs <- segs %>% rename(chr = Chromosome)
        if (!("startpos" %in% nm) && ("start" %in% nm)) segs <- segs %>% rename(startpos = start)
        if (!("endpos" %in% nm) && ("end" %in% nm)) segs <- segs %>% rename(endpos = end)
        # Ensure types
        segs <- segs %>% mutate(
            chr = toupper(gsub("^CHR", "", toupper(as.character(chr)))),
            startpos = as.numeric(startpos),
            endpos = as.numeric(endpos),
            nMajor = as.numeric(nMajor),
            nMinor = as.numeric(nMinor)
        )
    }

    # Ensure a sample column exists; fallback to provided sample_id or file name
    if (!("sample" %in% names(segs))) {
        default_sample <- if (!is.null(sample_id)) sample_id else tools::file_path_sans_ext(basename(segment_file))
        segs$sample <- default_sample
    }

    if (exclude_chrXY) {
        segs <- segs %>% filter(!chr %in% c("X", "Y"))
    }

    if (!is.null(sample_id) && ("sample" %in% names(segs))) {
        segs <- segs %>% filter(sample == sample_id)
    }

    segs <- segs %>%
        mutate(seg_size = endpos - startpos) %>%
        filter(seg_size >= min_seg_size)

    segs <- segs %>%
        mutate(chr = factor(chr, levels = chr_order)) %>%
        arrange(chr, startpos)

    # Collapse adjacent segments with identical copy numbers to eliminate visual gaps
    segs <- segs %>%
        group_by(chr, sample) %>%
        arrange(startpos) %>%
        mutate(
            # Create groups for consecutive segments with same copy numbers
            cn_group = cumsum(
                nMajor != lag(nMajor, default = dplyr::first(nMajor) + 1) |
                    nMinor != lag(nMinor, default = dplyr::first(nMinor) + 1)
            )
        ) %>%
        group_by(chr, sample, cn_group, nMajor, nMinor) %>%
        summarise(
            startpos = min(startpos),
            endpos = max(endpos),
            seg_size = endpos - startpos,
            .groups = "drop"
        ) %>%
        select(-cn_group) %>%
        arrange(chr, startpos)

    # Always offset alleles (create allele-specific long format with pair_id)
    segs_with_pairs <- segs %>%
        mutate(pair_id = row_number()) %>%
        pivot_longer(cols = c(nMajor, nMinor), names_to = "allele", values_to = "copy_number") %>%
        mutate(
            chr = factor(chr, levels = chr_order),
            offset_cn = case_when(
                allele == "nMajor" ~ copy_number + offset,
                allele == "nMinor" ~ copy_number - offset,
                TRUE ~ copy_number
            )
        ) %>%
        mutate(offset_cn = pmin(offset_cn, cn_cap))

    # Cap copy number at cn_cap (safety)
    segs_with_pairs <- segs_with_pairs %>% mutate(offset_cn = pmin(offset_cn, cn_cap))

    # Chromosome coordinates
    chr_lengths <- segs %>%
        group_by(chr) %>%
        summarise(chr_len = max(endpos), .groups = "drop") %>%
        mutate(chr = factor(chr, levels = chr_order)) %>%
        arrange(chr) %>%
        mutate(
            chr_len = as.numeric(chr_len),
            chr_start = c(0, head(cumsum(chr_len), -1))
        )

    # Join coordinates
    segs_with_pairs <- segs_with_pairs %>%
        left_join(chr_lengths, by = "chr") %>%
        mutate(
            start_genome = startpos + chr_start,
            end_genome = endpos + chr_start + 1 # Extend by 1 bp to eliminate visual gaps
        )

    if (any(is.na(segs_with_pairs$chr_start))) {
        stop("NA values in 'chr_start'. Check chromosome naming.")
    }

    # Axis labels and vertical lines
    chr_midpoints <- chr_lengths %>%
        mutate(mid = chr_start + chr_len / 2)

    chr_boundaries <- chr_lengths %>%
        mutate(boundary = chr_start + chr_len) %>%
        pull(boundary)

    # Build custom y-axis with "cn_cap+" label
    y_breaks <- seq(0, cn_cap, by = 1)
    y_labels <- as.character(y_breaks)
    y_labels[length(y_labels)] <- paste0(cn_cap, "+")

    # Plot
    ggplot(segs_with_pairs) +
        geom_segment(aes(
            x = start_genome, xend = end_genome,
            y = offset_cn, yend = offset_cn,
            color = allele
        ), linewidth = line_width) +
        geom_vline(xintercept = chr_boundaries, color = "gray90", linetype = "dashed", linewidth = 0.3) +
        scale_x_continuous(name = "Chromosome", breaks = chr_midpoints$mid, labels = chr_midpoints$chr) +
        scale_y_continuous(
            name = paste0("Allele-specific\nCopy Number"),
            breaks = y_breaks, labels = y_labels, expand = c(0.01, 0), limits = c(-0.2, cn_cap),
            minor_breaks = NULL
        ) +
        scale_color_manual(values = c(nMajor = nmaj_color, nMinor = nmin_color), name = "Allele") +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1)
        )
}

# CNA Heatmap
# Requires following function to create heatmap matrix from segment files
format_cna <- function(segments_raw, exclude_sex_chromosomes = TRUE, min_segment_size = 500000) {
    # Load data
    segs <- fread(segments_raw) %>% mutate(
        nMajor = as.numeric(nMajor),
        nMinor = as.numeric(nMinor)
    )

    # Create CN matrix for heatmap
    cna_seg <- segs
    colnames(cna_seg)[colnames(cna_seg) == "nMajor"] <- paste0(unique(cna_seg$sample), "_CN_major")
    colnames(cna_seg)[colnames(cna_seg) == "nMinor"] <- paste0(unique(cna_seg$sample), "_CN_minor")
    colnames(cna_seg)[colnames(cna_seg) == "startpos"] <- "start"
    colnames(cna_seg)[colnames(cna_seg) == "endpos"] <- "end"
    cna_seg <- cna_seg %>% dplyr::select(chr, start, end, paste0(unique(cna_seg$sample), "_CN_major"), paste0(unique(cna_seg$sample), "_CN_minor"))

    # Calculate percentage of genome affected by SCNA
    # Normalize chromosome names (in case of "chrX", "chrY", etc.)
    segs$chr <- gsub("^chr", "", segs$chr)
    # Convert chromosome column to character for filtering
    segs$chr <- as.character(segs$chr)
    # Optionally exclude sex chromosomes
    if (exclude_sex_chromosomes) {
        segs <- segs[!(segs$chr %in% c("X", "Y")), ]
    }
    # Add segment length
    segs$length <- segs$endpos - segs$startpos + 1
    # Filter out short segs
    segs <- segs[segs$length >= min_segment_size, ]
    # Calculate total copy number
    segs$total_cn <- segs$nMajor + segs$nMinor
    # Define altered segs (any deviation from diploid)
    altered_segs <- segs[segs$total_cn != 2, ]
    # Compute total and altered genome length
    total_length <- sum(segs$length)
    altered_length <- sum(altered_segs$length)
    # Compute SCNA burden percentage
    scna_percent <- (altered_length / total_length) * 100
    # Make it a df
    scna_percent_df <- data.frame(sample = unique(segs$sample), scna_percent = round(scna_percent, 2))

    # Return outputs as list
    return(list(scna_mat = cna_seg, scna_tbl = scna_percent_df))
}

# Function to plot a heatmap of copy number alterations across samples and chromosomes
cna_heatmap <- function(
    cna_input_csv,
    tmb_annot = TRUE, # whether to include TMB annotation
    order_by = "input" # "input", "tmb" or "scna"
    ) {
    # load libraries
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(ComplexHeatmap))
    suppressPackageStartupMessages(library(circlize))
    suppressPackageStartupMessages(library(GenomicRanges))

    # load input csv
    cna_input <- fread(cna_input_csv)

    # Replace NA in tumor_type with "Unknown"
    cna_input$tumor_type[is.na(cna_input$tumor_type)] <- "Unknown"

    #  extract sample names
    samples <- cna_input$sample

    # named list of segment files
    named_files <- setNames(as.list(cna_input$segments), samples)

    # Process cna segment data for heatmap
    # Apply to all patients
    list_cna_mat_percent <- purrr::map(named_files, format_cna)
    # Get list of cna matrices
    list_cna <- list_cna_mat_percent %>% purrr::map("scna_mat")

    # Make table of cna percent for all patients
    percent <- list_cna_mat_percent %>%
        purrr::map("scna_tbl") %>%
        bind_rows()
    # write table cna percent
    write.table(percent, file = "cna_percent.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

    # Get chromosome ranges from circlize package (read.chromInfo()$df) and make a GRanges object of the chromosome lengths
    chr_df <- read.chromInfo()$df
    chr_df <- chr_df[chr_df$chr %in% paste0("chr", c(1:22, "X")), ]
    chr_df$chr <- gsub("chr", "", chr_df$chr) # Only if needed to remove the "chr" prefix, depends on reference genome used
    chr_gr <- GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
    chr_gr
    # Create genomic windows of 1Mb across the genome
    chr_window <- EnrichedHeatmap::makeWindows(chr_gr, w = 1e6)
    chr_window
    chr_window

    # Convert CNA data to GRange objects matrix
    average_in_window <- function(window, gr, v, method = "weighted", empty_v = NA) {
        if (missing(v)) v <- rep(1, length(gr))
        if (is.null(v)) v <- rep(1, length(gr))
        if (is.atomic(v) && is.vector(v)) v <- cbind(v)

        v <- as.matrix(v)
        if (is.character(v) && ncol(v) > 1) {
            stop("`v` can only be a character vector.")
        }

        if (length(empty_v) == 1) {
            empty_v <- rep(empty_v, ncol(v))
        }

        u <- matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))

        mtch <- as.matrix(findOverlaps(window, gr))
        intersect <- pintersect(window[mtch[, 1]], gr[mtch[, 2]])
        w <- width(intersect)
        v <- v[mtch[, 2], , drop = FALSE]
        n <- nrow(v)

        ind_list <- split(seq_len(n), mtch[, 1])
        window_index <- as.numeric(names(ind_list))
        window_w <- width(window)

        if (is.character(v)) {
            for (i in seq_along(ind_list)) {
                ind <- ind_list[[i]]
                if (is.function(method)) {
                    u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
                } else {
                    tb <- tapply(w[ind], v[ind], sum)
                    u[window_index[i], ] <- names(tb[which.max(tb)])
                }
            }
        } else {
            if (method == "w0") {
                gr2 <- reduce(gr, min.gapwidth = 0)
                mtch2 <- as.matrix(findOverlaps(window, gr2))
                intersect2 <- pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])

                width_intersect <- tapply(width(intersect2), mtch2[, 1], sum)
                ind <- unique(mtch2[, 1])
                width_setdiff <- width(window[ind]) - width_intersect

                w2 <- width(window[ind])

                for (i in seq_along(ind_list)) {
                    ind <- ind_list[[i]]
                    x <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
                    u[window_index[i], ] <- (x * width_intersect[i] + empty_v * width_setdiff[i]) / w2[i]
                }
            } else if (method == "absolute") {
                for (i in seq_along(ind_list)) {
                    u[window_index[i], ] <- colMeans(v[ind_list[[i]], , drop = FALSE])
                }
            } else if (method == "weighted") {
                for (i in seq_along(ind_list)) {
                    ind <- ind_list[[i]]
                    u[window_index[i], ] <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
                }
            } else {
                if (is.function(method)) {
                    for (i in seq_along(ind_list)) {
                        ind <- ind_list[[i]]
                        u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
                    }
                } else {
                    stop("wrong method.")
                }
            }
        }

        return(u)
    }
    mat <- matrix(nrow = 3030, ncol = 0) # 3030 is the number of 1Mb windows in the human genome
    for (i in 1:length(list_cna)) {
        gr_cna <- GRanges(seqnames = list_cna[[i]]$chr, ranges = IRanges(list_cna[[i]]$start, list_cna[[i]]$end))
        v <- as.matrix(list_cna[[i]][, -(1:3)])
        num_mat <- average_in_window(chr_window, gr_cna, v, method = "absolute", empty_v = NA)
        mat <- cbind(mat, num_mat)
    }

    # Maintain chrom order in heatmap
    chr <- as.vector(seqnames(chr_window))
    chr_level <- paste0(c(1:22, "X"))
    chr <- factor(chr, levels = chr_level)

    # Process TMB info if tmb_annot is TRUE
    if (tmb_annot) {
        # convert missing NSM values to 0
        cna_input$nsm_snv[is.na(cna_input$nsm_snv)] <- 0
        cna_input$nsm_indel[is.na(cna_input$nsm_indel)] <- 0

        # Create matrix of SNV and INDEL counts
        snv_indel_mat <- cna_input %>%
            dplyr::select(nsm_snv, nsm_indel) %>%
            as.matrix()
        # Apply a cap of 200 to matrix (to sum SNV+INDEL but respect proportions)
        # Capping function
        cap_and_transform <- function(row, max_sum = 200) {
            if (all(is.na(row))) {
                return(row) # Return the row unchanged if it only contains NA
            }

            total <- sum(row, na.rm = TRUE) # Ignore NAs in the sum
            if (total > max_sum) {
                factor <- max_sum / total
                row[!is.na(row)] <- row[!is.na(row)] * factor # Scale only non-NA values
            }

            return(row)
        }
        # Apply capping function to each row
        capped_snv_indel_mat <- t(apply(snv_indel_mat, 1, cap_and_transform))

        # Create a vector for NSM numbers interleaving with empty strings, to annotate every other column (because each sample as two allele columns)
        cna_input$nsm <- cna_input$nsm_snv + cna_input$nsm_indel
        NSM_annot <- rep("", length(cna_input$nsm) * 2) # Preallocate a vector with empty strings
        NSM_annot[seq(2, length(NSM_annot), 2)] <- cna_input$nsm # Fill every odd index with the original values

        # Heatmap annot with TMB info
        column_ha <- HeatmapAnnotation(
            "TMB_txt" = anno_text(NSM_annot, location = 0.1, just = "left", gp = gpar(fontsize = 10, col = "grey30")),
            "TMB" = anno_barplot(capped_snv_indel_mat[rep(1:nrow(snv_indel_mat), each = 2), ], bar_width = 1, gp = gpar(fill = c("#00C5CD", "blue"), col = c("#00C5CD", "blue"))),
            "SCNA (% of genome)" = rep(percent$scna_percent, each = 2),
            "Tumor Group" = rep(cna_input$tumor_type, each = 2),
            annotation_name_side = "left", annotation_name_gp = gpar(col = "grey20", fontsize = 10),
            simple_anno_size = unit(0.5, "cm")
        )

        # SNV and INDEL legend
        TMB_legend <- Legend(labels = c("SNV", "INDEL"), title = "TMB", legend_gp = gpar(fill = c("#00C5CD", "blue")))
    } else {
        # Heatmap annot without TMB info
        column_ha <- HeatmapAnnotation(
            "SCNA (% of genome)" = rep(percent$scna_percent, each = 2),
            "Tumor Group" = rep(cna_input$tumor_type, each = 2),
            annotation_name_side = "left", annotation_name_gp = gpar(col = "grey20", fontsize = 10),
            simple_anno_size = unit(0.5, "cm")
        )
    }

    # cap values at 3 and round to nearest integer
    mat <- pmin(mat, 3)
    mat <- round(mat, 0)
    # Color function for heatmap
    heatmap_colors <- structure(c("#00868B", "grey98 ", "#B77FE2", "#7D26CD"), names = c("0", "1", "2", "3"))

    # Order samples by input, TMB or SCNA
    if (order_by == "input") {
        sample_order <- unique(cna_input$sample)
    } else if (order_by == "tmb") {
        sample_order <- cna_input %>%
            arrange(nsm) %>%
            distinct(sample)
    } else if (order_by == "scna") {
        sample_order <- percent %>%
            arrange(scna_percent) %>%
            distinct(sample)
    } else {
        stop("order_by must be either 'input', 'tmb' or 'scna'")
    }

    # Create heatmap
    hm <- Heatmap(mat,
        name = "Allelic Copy Number", top_annotation = column_ha, col = heatmap_colors,
        row_split = chr, cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(rep(unique(cna_input$sample), each = 2), levels = sample_order, labels = sample_order),
        show_column_names = FALSE,
        # column_title = NULL,
        column_title_rot = 90, column_title_gp = gpar(fontsize = 10),
        row_title_rot = 0, row_title_gp = gpar(fontsize = 7), border = TRUE,
        row_gap = unit(0, "points"), column_gap = unit(0, "points"), border_gp = gpar(col = "grey50", lwd = 0.5),
        heatmap_legend_param = list(at = c("0", "1", "2", "3"), labels = c("0", "1", "2", "3+"))
    )

    #  Draw heatmap (with or without TMB legend)
    if (tmb_annot) {
        draw(hm, annotation_legend_list = list(TMB_legend), merge_legend = FALSE)
    } else {
        draw(hm, merge_legend = FALSE)
    }
}
