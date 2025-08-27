#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
    cat("
How to use:
arg1: path to cfDNA bam
args2: path to normal bam
args3: patient name
args4: bed file
args5: path to ascat ref dir
args6: n cores to use
args7: sex ('XX' or 'XY')
args8: genome (hg19 or hg38)
args9: outdir
args10: skip normal (TRUE/FALSE)

      ")
    quit()
}

# Libraries
suppressPackageStartupMessages(library(ASCAT))

# Variables
tumourseqfile_path <- args[1]
normalseqfile_path <- args[2]
patient <- args[3]
bedfile_path <- args[4]
ascat_ref_path <- args[5]
threads <- args[6]
sex <- args[7] # 'XX' or 'XY'
genome <- args[8] #  hg19 or hg38
genome_number <- gsub("hg", "", genome)
bedfile_path <- gsub("hg\\d{2}", genome, bedfile_path)
outdir <- args[9]
skip_normal_process <- as.logical(args[10])  # converts "TRUE"/"FALSE" to logical

setwd(outdir)

# Tom Watkins' (tom.watkins@crick.ac.uk) function to run ASCAT in BAF-only mode for cfDNA samples (avoids over-segmentation due to more noisy sequencing)
# Requires an ascat.bc object that results from running the ascat.aspcf function:
# e.g.   ascat.bc = ascat.aspcf(ascat.bc)
# returns another ascat.bc object that has been modified so that only the BAF segmentation is applied to both BAF and LogR
# resulting copy number calls based only on the BAF segmentation when this ascat.bc object is used with ascat.runAscat function.
makeASCATObjectSegmentedOnlyOnBAFWithLogRCorrection <- function(ascat.bc) {
    print("Called makeASCATObjectSegmentedOnlyOnBAFWithLogRCorrection")
    # This is a version of the base R function rle that can also handle runs NAs which is not the default behaviour
    # Shouldn't be necessary for this task as no NAs should be present in the ascat object but I use it as it's more robust...
    rle_w_na_handling <- function(x) {
        if (!is.vector(x) && !is.list(x)) {
            stop("'x' must be an atomic vector")
        }
        n <- length(x)
        if (n == 0L) {
            return(structure(list(lengths = integer(), values = x),
                class = "rle"
            ))
        }
        IS_LOGIC <- FALSE
        #### BEGIN NEW SECTION PART 1 ####
        naRepFlag <- F
        if (any(is.na(x))) {
            naRepFlag <- T
            IS_LOGIC <- ifelse(typeof(x) == "logical", T, F)

            if (typeof(x) == "logical") {
                x <- as.integer(x)
                naMaskVal <- 2
            } else if (typeof(x) == "character") {
                naMaskVal <- paste(sample(c(letters, LETTERS, 0:9), 32, replace = T), collapse = "")
            } else {
                naMaskVal <- max(0, abs(x[!is.infinite(x)]), na.rm = T) + 1
            }

            x[which(is.na(x))] <- naMaskVal
        }
        #### END NEW SECTION PART 1 ####

        y <- x[-1L] != x[-n]
        i <- c(which(y), n)

        #### BEGIN NEW SECTION PART 2 ####
        if (naRepFlag) {
            x[which(x == naMaskVal)] <- NA
        }

        if (IS_LOGIC) {
            x <- as.logical(x)
        }
        #### END NEW SECTION PART 2 ####

        structure(list(lengths = diff(c(0L, i)), values = x[i]),
            class = "rle"
        )
    }

    tf.vec.seg.baf.snps <- rownames(ascat.bc[["Tumor_BAF"]]) %in% rownames(ascat.bc[["Tumor_BAF_segmented"]][[1]])
    # Adding the chromosome number to segmented baf here to make the runs take into account chromosome otherwise they
    # may span multiple chromosomes with the same segmented BAF value (e.g. 0.5)
    ascat.bc$SNPpos$Chromosome <- sub("X", 23, ascat.bc$SNPpos$Chromosome)
    ascat.bc$SNPpos$Chromosome <- sub("Y", 24, ascat.bc$SNPpos$Chromosome)
    # Here I take the
    rle.baf.pcf.df <- with(rle_w_na_handling(ascat.bc[["Tumor_BAF_segmented"]][[1]][, 1] + as.numeric(ascat.bc$SNPpos$Chromosome[tf.vec.seg.baf.snps])), data.frame(
        number = values,
        start = cumsum(lengths) - lengths + 1,
        end = cumsum(lengths)
    )[order(values), ])

    rle.baf.pcf.df <- rle.baf.pcf.df[with(rle.baf.pcf.df, order(start)), ]

    rle.baf.pcf.df$start_snp_name <- rownames(ascat.bc[["Tumor_BAF_segmented"]][[1]])[rle.baf.pcf.df$start]
    rle.baf.pcf.df$end_snp_name <- rownames(ascat.bc[["Tumor_BAF_segmented"]][[1]])[rle.baf.pcf.df$end]
    # Need to get the start end SNP indices as the tumor baf segmented vector only has a subset of the probes with logr
    # So need to match names to get equivalent runs in the logr.
    rle.baf.pcf.df$start_snp_idx <- which(rownames(ascat.bc[["Tumor_BAF"]]) %in% rle.baf.pcf.df$start_snp_name)
    rle.baf.pcf.df$end_snp_idx <- which(rownames(ascat.bc[["Tumor_BAF"]]) %in% rle.baf.pcf.df$end_snp_name)

    # Now we have the idxs of the starts and ends of the runs of segmented baf we can calculate the mean logr in the
    # of all SNPs in each run and assign that as the segmented logr value for that run of segmente BAF
    tum.logr.vec <- ascat.bc[["Tumor_LogR"]][, 1]
    tum.logr.segmented.vec <- ascat.bc[["Tumor_LogR_segmented"]][, 1]
    for (i in 1:nrow(rle.baf.pcf.df)) {
        mean.seg.logr <- mean(tum.logr.vec[rle.baf.pcf.df[i, ]$start_snp_idx:rle.baf.pcf.df[i, ]$end_snp_idx], na.rm = TRUE)
        # If there is a single value in this segment that isn't NA there should be a non-NA value
        # However, if there still is an NA value because there are no non-NA values in this segment it will be
        # set to zero...
        if (is.na(mean.seg.logr)) {
            print("Segment below did not contain a single value that was not NA, therefore setting the LogR to zero.")
            print(rle.baf.pcf.df[i, ])
            mean.seg.logr <- 0
        }
        tum.logr.segmented.vec[rle.baf.pcf.df[i, ]$start_snp_idx:rle.baf.pcf.df[i, ]$end_snp_idx] <- mean.seg.logr
    } # end for (i in 1:nrow(rle.baf.pcf.df)){

    # Bit of extra code I wrote just to check that the runs are now rougly comparable. As there are some logr values that
    # fall outside the BAF segmentation they will not be reset and there will always be more logr runs than BAF runs. The main
    # thing to make sure is that the runs from the segmented baf are now present in the logr segmented.
    # Commented out for now as only good for eyeballing interactively.
    # rle.logr.df <- with(rle_w_na_handling(tum.logr.segmented.vec), data.frame(number = values,
    #                                                                           start = cumsum(lengths) - lengths + 1,
    #                                                                           end = cumsum(lengths))[order(values),])
    # rle.logr.df <- rle.logr.df[with(rle.logr.df, order(start)),]

    # Set the values that cause the BAF segmentation to be applied to the LogR in the ascat object
    ascat.bc[["Tumor_LogR_segmented"]][, 1] <- tum.logr.segmented.vec

    rle.logr.pcf.df <- with(rle_w_na_handling(ascat.bc[["Tumor_LogR_segmented"]][, 1]), data.frame(
        number = values,
        start = cumsum(lengths) - lengths + 1,
        end = cumsum(lengths)
    )[order(values), ])
    rle.logr.pcf.df <- rle.logr.pcf.df[with(rle.logr.pcf.df, order(start)), ]

    rle.baf.pcf.df$start_end_name <- paste0(rle.baf.pcf.df$start_snp_idx, "_", rle.baf.pcf.df$end_snp_idx)
    rle.logr.pcf.df$start_end_name <- paste0(rle.logr.pcf.df$start, "_", rle.logr.pcf.df$end)
    rle.logr.pcf.df$baf_segment <- FALSE

    rle.logr.pcf.df[rle.logr.pcf.df$start_end_name %in% rle.baf.pcf.df$start_end_name, ]$baf_segment <- TRUE
    # Now we need to add in that the chromsome of these segments so that we don't take a value corresponding to a different one.
    rle.logr.pcf.df$chrom <- NA
    rle.logr.pcf.df$num_pos <- NA

    old.rle.logr.pcf.df <- rle.logr.pcf.df
    print(rle.logr.pcf.df)

    i <- 1
    curr.nrow.rle.logr.pcf.df <- nrow(rle.logr.pcf.df)

    while (i <= curr.nrow.rle.logr.pcf.df) {
        seg.present.on.chrom <- unique(as.numeric(ascat.bc$SNPpos$Chromosome[as.numeric(rle.logr.pcf.df[i, ]$start):as.numeric(rle.logr.pcf.df[i, ]$end)]))

        if (length(seg.present.on.chrom) > 1) {
            print(seg.present.on.chrom)
            print(rle.logr.pcf.df[i, ])

            if (length(seg.present.on.chrom) == 2) {
                # stop()
                # So to fix this segment we need to find the chromosome boundary that's crossed in terms of SNP positons
                # we know the lower number chromosome and higher number (the only issue is there migh tbe X and Y chroms here, but let's ignore for now)
                end.overlap.chrom <- seg.present.on.chrom[1]
                start.overlap.chrom <- seg.present.on.chrom[2]
                chr.seg.logr.df <- cbind(ascat.bc$SNPpos, ascat.bc$Tumor_LogR_segmented, 1:length(ascat.bc$Tumor_LogR_segmented), stringsAsFactors = FALSE)
                colnames(chr.seg.logr.df)[4] <- "snp_idx"
                # Ok let's look up the corresponding snp positions for both of these cases.
                chr.seg.logr.df$start_chr_idx <- FALSE
                chr.seg.logr.df[which(!duplicated(chr.seg.logr.df$chrom)), ]$start_chr_idx <- TRUE
                # get last index of each chromosome
                chr.seg.logr.df$end_chr_idx <- FALSE
                chr.seg.logr.df[length(chr.seg.logr.df$chrom) - match(unique(chr.seg.logr.df$chrom), rev(chr.seg.logr.df$chrom)) + 1, ]$end_chr_idx <- TRUE

                end.overlap.chrom.df <- chr.seg.logr.df[chr.seg.logr.df$chrom %in% end.overlap.chrom & chr.seg.logr.df$end_chr_idx %in% TRUE, , drop = FALSE]
                start.overlap.chrom.df <- chr.seg.logr.df[chr.seg.logr.df$chrom %in% start.overlap.chrom & chr.seg.logr.df$start_chr_idx %in% TRUE, , drop = FALSE]

                # So now we need to get this row out of the rle.logr.pc.df
                all.curr.idx.vec <- 1:nrow(rle.logr.pcf.df)
                new.idx.vec <- all.curr.idx.vec[!all.curr.idx.vec %in% i]
                # So now we have sufficient information to split the

                removed.pcf.row.df <- rle.logr.pcf.df[i, , drop = FALSE]

                new.first.pcf.row.df <- removed.pcf.row.df
                new.first.pcf.row.df$end <- end.overlap.chrom.df$snp_idx
                new.first.pcf.row.df$start_end_name <- paste0(new.first.pcf.row.df$start, "_", new.first.pcf.row.df$end)
                new.first.pcf.row.df$chrom <- end.overlap.chrom
                new.first.pcf.row.df$num_pos <- length(new.first.pcf.row.df$start:new.first.pcf.row.df$end)
                new.second.pcf.row.df <- removed.pcf.row.df
                new.first.pcf.row.df$chrom <- start.overlap.chrom
                new.second.pcf.row.df$start <- start.overlap.chrom.df$snp_idx
                new.second.pcf.row.df$num_pos <- length(new.second.pcf.row.df$start:new.second.pcf.row.df$end)
                new.second.pcf.row.df$start_end_name <- paste0(new.second.pcf.row.df$start, "_", new.second.pcf.row.df$end)
                # Now add these back in to the data frame..

                rle.logr.pcf.df <- rbind(new.first.pcf.row.df, new.second.pcf.row.df, rle.logr.pcf.df[new.idx.vec, ], stringsAsFactors = FALSE)
                rle.logr.pcf.df <- rle.logr.pcf.df[with(rle.logr.pcf.df, order(start, end)), ]
                # reset limit on loop
                curr.nrow.rle.logr.pcf.df <- nrow(rle.logr.pcf.df)
            } else {
                stop("There is a segment that spans at least three chromosomes...")
            }
        } else { # end if(length(seg.present.on.chrom) > 1) {

            rle.logr.pcf.df[i, ]$chrom <- seg.present.on.chrom
            rle.logr.pcf.df[i, ]$num_pos <- length(rle.logr.pcf.df[i, ]$start:rle.logr.pcf.df[i, ]$end)

            i <- i + 1
        }
    } # end for (i in 1:nrow(rle.logr.pcf.df)) {


    # This is the dataframe that we will use to test the result of our joining to attempt to make sure we haven't messed up the segmentation...
    orig.rle.logr.pcf.df <- rle.logr.pcf.df
    # Ok so now we need to go through these dataframes and make sure that the segmented logr is incorporated into an expaned BAF segment:
    # So now we need to flatten these
    j <- 1
    removed.row.df <- NA
    while (j <= nrow(rle.logr.pcf.df)) {
        # Taken the decision to join non-baf segments to their nearest neighbour, even if they are also a non-baf segment.
        # Choosing to join the segment to its neighbour with the closest logr value
        if (rle.logr.pcf.df[j, ]$baf_segment %in% FALSE) {
            # This will hold the mean logr values of the preceding and following segments on the same chromosomse or be NA if they don't exit.
            logr.comparison.vec <- c(NA, NA) # first entry is for preceding, second is for following segment.
            if ((j - 1) != 0 && rle.logr.pcf.df[(j - 1), ]$chrom == rle.logr.pcf.df[j, ]$chrom) {
                logr.comparison.vec[1] <- rle.logr.pcf.df[(j - 1), ]$number
            } # end if ( (j - 1) != 0  && rle.logr.pcf.df[ (j - 1) ,]$chrom == rle.logr.pcf.df[j,]$chrom) {

            if ((j + 1) <= nrow(rle.logr.pcf.df) && rle.logr.pcf.df[(j + 1), ]$chrom == rle.logr.pcf.df[j, ]$chrom) {
                logr.comparison.vec[2] <- rle.logr.pcf.df[(j + 1), ]$number
            } # end if ( (j - 1) != 0  && rle.logr.pcf.df[ (j - 1) ,]$chrom == rle.logr.pcf.df[j,]$chrom) {
            # which is closest to the current seg's logr value?
            logr.difference.vec <- abs(logr.comparison.vec - rle.logr.pcf.df[j, ]$number)

            # Was there at least one applicable adjacent segment to which this segment may be joined?
            if (!all(is.na(logr.difference.vec))) {
                # Should we join this segment to the preceding or following segment?
                prec.or.follow.seg <- which.min(logr.difference.vec)
                # This row is going to be joined, so we're going to have to remove it from the datframe... let's get the indices
                curr.rows.idx.vec <- 1:nrow(rle.logr.pcf.df)
                new.rows.idx.vec <- curr.rows.idx.vec[!curr.rows.idx.vec %in% j]
                removed.row.df <- rle.logr.pcf.df[j, , drop = FALSE]
                rle.logr.pcf.df <- rle.logr.pcf.df[new.rows.idx.vec, ]

                # Because we delete the current row, we keep the baf_segment classification of the new.
                if (prec.or.follow.seg == 1) {
                    # join the with the preceding segment
                    # So let's save the relevant info from the segment we're going to delete (ie the "middle one")
                    rle.logr.pcf.df[(j - 1), ]$end <- removed.row.df[1, ]$end
                    rle.logr.pcf.df[(j - 1), ]$start_end_name <- paste0(rle.logr.pcf.df[(j - 1), ]$start, "_", rle.logr.pcf.df[(j - 1), ]$end)
                    rle.logr.pcf.df[(j - 1), ]$num_pos <- length(as.numeric(rle.logr.pcf.df[(j - 1), ]$start):as.numeric(rle.logr.pcf.df[(j - 1), ]$end))
                    # So now we need to update the mean logr value for these segments...

                    mean.seg.logr <- mean(tum.logr.vec[rle.logr.pcf.df[(j - 1), ]$start:rle.logr.pcf.df[(j - 1), ]$end], na.rm = TRUE)
                    # If there is a single value in this segment that isn't NA there should be a non-NA value
                    # However, if there still is an NA value because there are no non-NA values in this segment it will be
                    # set to zero...
                    if (is.na(mean.seg.logr)) {
                        print("Segment below did not contain a single value that was not NA, therefore setting the LogR to zero.")
                        print(rle.logr.pcf.df[j, ])
                        mean.seg.logr <- 0
                    }
                    tum.logr.segmented.vec[rle.logr.pcf.df[(j - 1), ]$start:rle.logr.pcf.df[(j - 1), ]$end] <- mean.seg.logr
                    # So now in this case we need to revisit this segment again in case it to needs to be joined?
                    # No, it cannot be a logR segment as it would have been joined in that case, so leave j the same
                    # and the next segment will now be indexed by j and examined.
                } else if (prec.or.follow.seg == 2) { # end if (prec.or.follow.seg == 1) {
                    # join with the following segment
                    # So let's save the relevant info from the segment we're going to delete (ie the "middle one")
                    # So since we're removing our current row, the next idx in the original dataframe will
                    # become the row that we're working on.
                    # Start should be updated...
                    rle.logr.pcf.df[j, ]$start <- removed.row.df[1, ]$start
                    rle.logr.pcf.df[j, ]$start_end_name <- paste0(rle.logr.pcf.df[j, ]$start, "_", rle.logr.pcf.df[j, ]$end)
                    rle.logr.pcf.df[j, ]$num_pos <- length(as.numeric(rle.logr.pcf.df[j, ]$start):as.numeric(rle.logr.pcf.df[j, ]$end))
                    # So now we need to update the mean logr value for these segments...

                    mean.seg.logr <- mean(tum.logr.vec[rle.logr.pcf.df[j, ]$start:rle.logr.pcf.df[j, ]$end], na.rm = TRUE)
                    # If there is a single value in this segment that isn't NA there should be a non-NA value
                    # However, if there still is an NA value because there are no non-NA values in this segment it will be
                    # set to zero...
                    if (is.na(mean.seg.logr)) {
                        print("Segment below did not contain a single value that was not NA, therefore setting the LogR to zero.")
                        print(rle.logr.pcf.df[j, ])
                        mean.seg.logr <- 0
                    } # end if (is.na(mean.seg.logr)) {
                    tum.logr.segmented.vec[rle.logr.pcf.df[j, ]$start:rle.logr.pcf.df[j, ]$end] <- mean.seg.logr

                    # Now if this segment is a logr segment it will need to be revisited so let's leave j the same
                    # again.
                } # end } else if (prec.or.follow.seg == 2)
            } else { # end  if ( ! all(is.na(logr.difference.vec)) ) {
                # No segment that this segment could be joined to, increment j and move on
                j <- j + 1
            } # end else {#end  if ( ! all(is.na(logr.difference.vec)) ) {
        } else { # end  if (rle.logr.pcf.df[j,]$baf_segment %in% FALSE) {
            # The current segment is a BAF segment and should not be merged, increment j and move on
            j <- j + 1
        } # end  else {#end  if (rle.logr.pcf.df[j,]$baf_segment %in% FALSE) {

        # A good test to see if this joining worked (or one of them) is to check that the same number of SNP positions are represented.
        if (sum(orig.rle.logr.pcf.df$num_pos) != sum(rle.logr.pcf.df$num_pos)) {
            print("old:")
            print(orig.rle.logr.pcf.df)
            print("new:")
            print(rle.logr.pcf.df)
            print("segment to merge:")
            print(removed.row.df)
            print("current index:")
            print(j)
            print("Total old:")
            print(sum(orig.rle.logr.pcf.df$num_pos))
            print("Total new:")
            print(sum(rle.logr.pcf.df$num_pos))
            print("ERROR: Different total number of SNPs in the BAF segmentation after joining remaining logr only segments...")
            stop()
        } else { # end if ( sum(orig.rle.logr.pcf.df$num_pos) != sum(rle.logr.pcf.df$num_pos) ){

            print("Successfully merged segment:")
            print(removed.row.df)
            print("current index:")
            print(j)
            print("Total old SNPs:")
            print(sum(orig.rle.logr.pcf.df$num_pos))
            print("Total new SNPs:")
            print(sum(rle.logr.pcf.df$num_pos))
        } # end else {#end if ( sum(orig.rle.logr.pcf.df$num_pos) != sum(rle.logr.pcf.df$num_pos) ){
    } # end while (j <= (rle.logr.pcf.df) ) {

    # A good test to see if this joining worked (or one of them) is to check that the same number of SNP positions are represented.
    if (sum(orig.rle.logr.pcf.df$num_pos) != sum(rle.logr.pcf.df$num_pos)) {
        print(sum(orig.rle.logr.pcf.df$num_pos))
        print(sum(rle.logr.pcf.df$num_pos))
        print("ERROR: Different total number of SNPs in the BAF segmentation after joining remaining logr only segments...")
        # stop()
    } # if ( sum(orig.rle.logr.pcf.df$num_pos) != sum(rle.logr.pcf.df$num_pos) ){
    # Also test that all the chromosomes present in the original are present
    if (length(unique(orig.rle.logr.pcf.df$chrom)) != length(unique(rle.logr.pcf.df$chrom))) {
        print("ERROR: Different total number of chromosomes in the BAF segmentation after joining remaining logr only segments...")
        # stop()
    } # if ( length(unique(orig.rle.logr.pcf.df$chrom)) != length(unique(rle.logr.pcf.df$chrom)) ){
    # Are the same number of BAF segments present as in the original dataframe?
    if (nrow(orig.rle.logr.pcf.df[orig.rle.logr.pcf.df$baf_segment %in% TRUE, ]) != nrow(rle.logr.pcf.df[rle.logr.pcf.df$baf_segment %in% TRUE, ])) {
        print("ERROR: Different number of BAF segments (defined by the initial segmented BAF) after joining remaining logr only segments...")
        # stop()
    } # end if ( nrow(orig.rle.logr.pcf.df[orig.rle.logr.pcf.df$baf_segment %in% TRUE,]) != nrow(rle.logr.pcf.df[rle.logr.pcf.df$baf_segment %in% TRUE,]) ){

    # if ( any(rle.logr.pcf.df$) ){

    #   print("ERROR: Different number of BAF segments (defined by the initial segmented BAF) after joining remaining logr only segments...")
    #   #stop()

    # }
    # Test whether all the start and end snp positions are consecutive in all the segments.
    # Test whether start and end of each chromosome are the same as in the original segmentation:
    # maybe imlement later...
    # If these tests are successful then we can update the logr vector in the actual ascat.bc object...
    # There is also the potenital issue of segments that span chromosome arms... that I haven't addressed at all. Can revisit this.
    ascat.bc[["Tumor_LogR_segmented"]][, 1] <- tum.logr.segmented.vec

    ascat.bc$SNPpos$Chromosome <- sub("23", "X", ascat.bc$SNPpos$Chromosome)
    ascat.bc$SNPpos$Chromosome <- sub("24", "Y", ascat.bc$SNPpos$Chromosome)

    return(ascat.bc)
} # end makeASCATObjectSegmentedOnlyOnBAF <- function(ascat.bc) {


# References
alleles <- paste0(ascat_ref_path, "/G1000_allelesAll_hg", genome_number, "/G1000_alleles_hg", genome_number, "_chr")
loci <- paste0(ascat_ref_path, "/G1000_lociAll_hg", genome_number, "/G1000_loci_hg", genome_number, "_chr") # Need to change chromosome name!! for i in {1..22} X; do sed -i 's/^/chr/' G1000_loci_hg", genome_number, "_chr${i}.txt; done
gc <- paste0(ascat_ref_path, "/GC_G1000_hg", genome_number, ".txt")
rt <- paste0(ascat_ref_path, "/RT_G1000_hg", genome_number, ".txt")

ascat.prepareHTS(
    tumourseqfile = tumourseqfile_path,
    normalseqfile = normalseqfile_path,
    tumourname = paste0(patient, "_cfDNA"),
    normalname = paste0(patient, "_Normal"),
    allelecounter_exe = "/opt/bin/alleleCounter",
    alleles.prefix = alleles,
    loci.prefix = loci,
    gender = sex,
    genomeVersion = genome,
    nthreads = threads,
    BED_file = bedfile_path,
    chrom_names = c(1:22, "X"),
	skip_allele_counting_tumour = FALSE,
	skip_allele_counting_normal = skip_normal_process) # set to TRUE if already done

ascat.bc <- ascat.loadData(Tumor_LogR_file = paste0(patient, "_cfDNA_tumourLogR.txt"), Tumor_BAF_file = paste0(patient, "_cfDNA_tumourBAF.txt"), Germline_LogR_file = paste0(patient, "_cfDNA_normalLogR.txt"), Germline_BAF_file = paste0(patient, "_cfDNA_normalBAF.txt"), gender = sex, genomeVersion = genome)
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = gc, replictimingfile = rt)
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc <- ascat.aspcf(ascat.bc)
ascat.bc <- makeASCATObjectSegmentedOnlyOnBAFWithLogRCorrection(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output <- ascat.runAscat(ascat.bc, gamma = 1, write_segments = T)
QC <- ascat.metrics(ascat.bc, ascat.output)
save(ascat.bc, ascat.output, QC, file = "ASCAT_objects_cfDNA.Rdata")
