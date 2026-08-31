# Calculate strand bias from paired watson and crick polymerase usage wig files.
#
# Strand bias is calculated for each genomic bin as:
#
#   (Watson - Crick) / (Watson + Crick)
#
# The resulting values are exported as variable-step WIG files.


library(dplyr)
library(rtracklayer)

# Define the list of input files and output directory

watson_list <- list(
  "./data/wig/pol-h-f18a_pol-usage_watson_w1000_MA30.bw")

crick_list <- list(
  "./data/wig/pol-h-f18a_pol-usage_crick_w1000_MA30.bw")

output.dir <- "./"

# Process the data

for (i in seq_along(watson_list)) {
  f.file <- watson_list[[i]]
  r.file <- crick_list[[i]]
  
  # Establish parameters and check for errors
  
  pol = sub(".*\\/(.*)\\.a1.*", "\\1", f.file)
  r.pol = sub(".*\\/(.*)\\.a1.*", "\\1", r.file)
  
  if (pol != r.pol) {
    stop("Error: polymerase must be the same.")
  }
  
  bin.width = sub(".*_w([0-9]+)_MA.*", "\\1", f.file)
  r.bin.width = sub(".*_w([0-9]+)_MA.*", "\\1", r.file)
  
  if (bin.width != r.bin.width) {
    stop("Error: bin width must be the same value.")
  }
  
  MA = sub(".*_MA([0-9]+)_.*", "\\1", f.file)
  r.MA = sub(".*_MA([0-9]+)_.*", "\\1", r.file)
  
  if (MA != r.MA) {
    stop("Error: MA must be the same value.")
  }
  
  if(!dir.exists(output.dir))dir.create(output.dir)
  
  if (output.dir != "") {
    file.name <- paste(output.dir, "/", pol, ".w", bin.width, "_MA", MA, ".str-bias.wig", sep = "")
  } else {
    file.name <- paste(pol, ".w", bin.width, "_MA", MA, ".str-bias.wig", sep = "")
  }
  
  track.name <- paste(pol, "_MA", MA, ".str-bias.wig", sep = "")
  current_date <- Sys.Date()
  track.decription <- paste("Analysis by L. Bainbridge; ", current_date)
  additional_params <- "visibility=full autoScale=off color=214,43,59 yLineOnOff=on yLineMark=0 viewLimits=-0.5:0.5 priority=10"
  header <- paste("track type=wiggle_0 name=\"", track.name, "\" description=\"", track.decription, "\"", " ", additional_params, sep="")
  
  
  #Read count files
  column_names = c("chro", "pos", "score")
  
  f.wig_file <- import(f.file)
  f.wig_file <- as.data.frame(f.wig_file)
  f.wig_file <- f.wig_file[, c(1, 2, 6)]
  colnames(f.wig_file) <- column_names
  
  r.wig_file <- import(r.file)
  r.wig_file <- as.data.frame(r.wig_file)
  r.wig_file <- r.wig_file[, c(1, 2, 6)]
  colnames(r.wig_file) <- column_names
  
  cat("wig files loaded successfully\n")
  
  # Merge polymerase WIG files
  merged_df <- merge(f.wig_file, r.wig_file, by = c("chro", "pos"), all = TRUE) %>%
    arrange(chro, pos)
  
  merged_df <- merged_df[complete.cases(merged_df), ]
  
  
  # Calculate strand bias
  merged_df$SB <- (merged_df$score.x - merged_df$score.y)/(merged_df$score.x + merged_df$score.y)
  
  
  # Export wig file
  cat("exporting polymerase track wig file...\n")
  
  writeLines(header, file.name)
  for (chromosome in unique(merged_df$chro)) {
    cat(paste("variableStep chrom=", chromosome, " span=1\n", sep=""), file = file.name, append = TRUE)
    subset_df <- merged_df[merged_df$chro == chromosome, ]
    for (i in 1:nrow(subset_df)) {
      cat(paste(subset_df$pos[i], subset_df$SB[i], sep="\t"), "\n", file = file.name, append = TRUE)
    }
    cat(chromosome, "processing complete\n")
  }
  
  cat("export complete\n")
  
}