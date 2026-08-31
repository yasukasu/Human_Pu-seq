# Generate a genome-wide gene-expression track from a BED file containing
# gene coordinates and FPKM values.
#
# For each 1-kb genomic window, the script sums the FPKM values of all
# overlapping genes. The resulting signal is smoothed using a centred
# moving average extending 30 windows (30 kb) in each direction.
#
# Note: this represents the expression of genes spanning each position,
# not RNA-seq read coverage.

library(data.table)
library(GenomicRanges)

# Input, output and analysis parameters

bed <- fread("./data/bed/hct116.FPKM.TSS-TTS.gold.bed", 
             col.names = c("chr", "start", "end", "id", "fpkm"))

chrom_sizes <- fread("./data/genome/GRCh38.chrom.sizes", 
                     col.names = c("chr", "size"))

output_name <- "gene_fpkm_signal_smoothed_MA30.wig"
output_dir <- "."

# Normalize BED coordinates
bed[, `:=`(start = pmin(start, end), end = pmax(start, end))]

# Filter for FPKM > 0
bed <- bed[fpkm > 0]

# Convert to GRanges
gr_genes <- GRanges(seqnames = bed$chr,
                    ranges = IRanges(start = bed$start + 1, end = bed$end),  # BED is 0-based
                    score = bed$fpkm)

# Build 1 kb genome-wide windows
window_size <- 1000
gr_windows <- GRangesList(lapply(1:nrow(chrom_sizes), function(i) {
  chr <- chrom_sizes$chr[i]
  seqlen <- chrom_sizes$size[i]
  GRanges(seqnames = chr,
          ranges = IRanges(start = seq(1, seqlen, by = window_size),
                           width = window_size))
})) |> unlist()

# Find overlaps
hits <- findOverlaps(gr_windows, gr_genes)

# Aggregate FPKM values per window
score <- numeric(length(gr_windows))
score_df <- data.table(query = queryHits(hits), fpkm = mcols(gr_genes)$score[subjectHits(hits)])
score_sum <- score_df[, .(fpkm_sum = sum(fpkm)), by = query]
score[score_sum$query] <- score_sum$fpkm_sum


# Smooth the data using a centred moving average (±30 bins)

span <- 61
half_span <- floor(span / 2)
kernel <- rep(1 / span, span)
smoothed_score <- numeric(length(score))

for (chr in unique(as.character(seqnames(gr_windows)))) {
  idx <- which(as.character(seqnames(gr_windows)) == chr)
  chr_scores <- score[idx]
  n <- length(chr_scores)
  
  if (n >= span) {
    # Calculate the complete 61-bin moving average.
    chr_smoothed <- as.numeric(
      stats::filter(chr_scores, kernel, sides = 2)
    )
    
    # At chromosome boundaries, average only the available bins.
    edge_positions <- which(is.na(chr_smoothed))
    
    for (i in edge_positions) {
      left <- max(1, i - half_span)
      right <- min(n, i + half_span)
      
      chr_smoothed[i] <- mean(chr_scores[left:right])
    }
  } else {
    # For short chromosomes, average all available neighbouring bins.
    chr_smoothed <- vapply(
      seq_len(n),
      function(i) {
        left <- max(1, i - half_span)
        right <- min(n, i + half_span)
        
        mean(chr_scores[left:right])
      },
      numeric(1)
    )
  }
  
  smoothed_score[idx] <- chr_smoothed
}

# Write WIG file
wig_smooth <- character()
for (chr in unique(as.character(seqnames(gr_windows)))) {
  idx <- which(as.character(seqnames(gr_windows)) == chr)
  wig_smooth <- c(wig_smooth,
                  sprintf("fixedStep chrom=%s start=1 step=%d span=%d", chr, window_size, window_size),
                  format(smoothed_score[idx], scientific = FALSE, trim = TRUE)
  )
}

output_file = file.path(output_dir, output_name)
if(!dir.exists(output.dir))dir.create(output.dir, recursive = TRUE)

writeLines(wig_smooth, con = output_file)

