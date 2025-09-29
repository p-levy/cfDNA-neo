#!/usr/bin/env Rscript

# Function to convert PURPLE segment TSV format to ASCAT segments_raw.txt format
# 
# This function takes a PURPLE segment file or data frame and converts it to the ASCAT format
# Expected PURPLE columns: chromosome, start, end, majorAlleleCopyNumber, minorAlleleCopyNumber
# Output ASCAT columns: sample, chr, startpos, endpos, nMajor, nMinor, nAraw, nBraw
#
# @param purple_input Path to PURPLE segment TSV file OR data frame with PURPLE data
# @param output_file Path for output ASCAT segments_raw.txt file (optional, if NULL returns data frame only)
# @param sample_name Sample name to use in output (if NULL, uses input filename or "sample")
# @param min_seg_size Minimum segment size to keep (default: 1e6 bp)
# @param exclude_chrXY Whether to exclude sex chromosomes (default: FALSE)
#
# @return Returns the converted data frame (and optionally writes to file)
convert_purple_to_ascat <- function(
    purple_input, 
    output_file = NULL, 
    sample_name = NULL,
    min_seg_size = 1e6,
    exclude_chrXY = FALSE) {
    
    # Load required libraries
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(dplyr))
    
    # Read PURPLE data (either from file or data frame)
    if (is.character(purple_input)) {
        cat("Reading PURPLE file:", purple_input, "\n")
        purple_data <- fread(purple_input)
        input_filename <- purple_input
    } else if (is.data.frame(purple_input)) {
        cat("Using provided PURPLE data frame\n")
        purple_data <- as.data.table(purple_input)
        input_filename <- "purple_data"
    } else {
        stop("purple_input must be either a file path (character) or a data frame")
    }
    
    # Check required columns
    required_cols <- c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber")
    missing_cols <- setdiff(required_cols, names(purple_data))
    if (length(missing_cols) > 0) {
        # Try alternative column names
        alt_names <- list(
            chromosome = c("Chromosome", "chr"),
            start = c("Start", "startpos", "minStart"),
            end = c("End", "endpos", "maxStart"),
            majorAlleleCopyNumber = c("nMajor", "major_cn"),
            minorAlleleCopyNumber = c("nMinor", "minor_cn")
        )
        
        for (col in missing_cols) {
            found <- FALSE
            for (alt in alt_names[[col]]) {
                if (alt %in% names(purple_data)) {
                    names(purple_data)[names(purple_data) == alt] <- col
                    found <- TRUE
                    break
                }
            }
            if (!found) {
                stop("Missing required column: ", col, " (and alternatives: ", 
                     paste(alt_names[[col]], collapse = ", "), ")")
            }
        }
    }
    
    # Determine sample name
    if (is.null(sample_name)) {
        if (is.character(purple_input)) {
            sample_name <- tools::file_path_sans_ext(basename(purple_input))
            sample_name <- gsub("\\.purple\\.segment$", "", sample_name)
        } else {
            sample_name <- "sample"
        }
    }
    
    cat("Converting to ASCAT format for sample:", sample_name, "\n")
    
    # Convert to ASCAT format
    ascat_data <- purple_data %>%
        mutate(
            sample = sample_name,
            chr = as.character(gsub("^chr", "", chromosome)),  # Remove "chr" prefix if present
            startpos = as.numeric(start),
            endpos = as.numeric(end),
            nMajor = round(as.numeric(majorAlleleCopyNumber)),
            nMinor = round(as.numeric(minorAlleleCopyNumber)),
            # ASCAT raw values (nAraw, nBraw) - use rounded values as approximation
            nAraw = as.numeric(majorAlleleCopyNumber),
            nBraw = as.numeric(minorAlleleCopyNumber)
        ) %>%
        select(sample, chr, startpos, endpos, nMajor, nMinor, nAraw, nBraw)
    
    # Filter by segment size
    if (min_seg_size > 0) {
        initial_count <- nrow(ascat_data)
        ascat_data <- ascat_data %>%
            mutate(seg_size = endpos - startpos) %>%
            filter(seg_size >= min_seg_size) %>%
            select(-seg_size)
        cat("Filtered segments by size (>=", min_seg_size, "bp):", 
            initial_count, "->", nrow(ascat_data), "\n")
    }
    
    # Exclude sex chromosomes if requested
    if (exclude_chrXY) {
        initial_count <- nrow(ascat_data)
        ascat_data <- ascat_data %>% filter(!chr %in% c("X", "Y"))
        cat("Excluded sex chromosomes:", initial_count, "->", nrow(ascat_data), "\n")
    }
    
    # Sort by chromosome and position
    chr_order <- c(as.character(1:22), "X", "Y")
    ascat_data <- ascat_data %>%
        mutate(chr = factor(chr, levels = chr_order)) %>%
        arrange(chr, startpos) %>%
        mutate(chr = as.character(chr))
    
    # Write output file if specified
    if (!is.null(output_file)) {
        cat("Writing ASCAT file:", output_file, "\n")
        write.table(ascat_data, file = output_file, sep = "\t", row.names = FALSE, 
                    quote = FALSE, col.names = TRUE)
    }
    
    cat("Conversion complete! Converted", nrow(ascat_data), "segments\n")
    cat("Sample:", sample_name, "\n")
    cat("Chromosomes:", paste(unique(ascat_data$chr), collapse = ", "), "\n")
    
    # Return data
    return(ascat_data)
}

# Example usage function
example_usage <- function() {
    cat("Example usage:\n")
    cat("# Convert single file\n")
    cat("convert_purple_to_ascat(\n")
    cat("  purple_file = 'sample.purple.segment.tsv',\n")
    cat("  output_file = 'sample.segments_raw.txt',\n")
    cat("  sample_name = 'sample_01'\n")
    cat(")\n\n")
    
    cat("# Convert with filtering\n")
    cat("convert_purple_to_ascat(\n")
    cat("  purple_file = 'sample.purple.segment.tsv',\n")
    cat("  output_file = 'sample.segments_raw.txt',\n")
    cat("  min_seg_size = 500000,  # 500kb minimum\n")
    cat("  exclude_chrXY = TRUE\n")
    cat(")\n")
}

# Print usage if script is run directly
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) == 0) {
        cat("PURPLE to ASCAT segment converter\n")
        cat("================================\n\n")
        example_usage()
        cat("\nCommand line usage:\n")
        cat("Rscript convert_purple_to_ascat.R input.purple.segment.tsv output.segments_raw.txt [sample_name]\n")
    } else if (length(args) >= 2) {
        sample_name <- if (length(args) >= 3) args[3] else NULL
        convert_purple_to_ascat(args[1], args[2], sample_name)
    } else {
        cat("Error: Need at least input and output file paths\n")
        cat("Usage: Rscript convert_purple_to_ascat.R input.purple.segment.tsv output.segments_raw.txt [sample_name]\n")
    }
}