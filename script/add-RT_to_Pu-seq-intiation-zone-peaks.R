# Add replication-timing values to initiation-site coordinates
#
# For each initiation site, this script assigns the replication-timing (RT)
# value from the nearest RT measurement on the same chromosome. The resulting
# bedGraph contains 1-bp initiation zone peaks, with RT used as the score.

library(rtracklayer)
library(dplyr)

# Define input and output parameters

bedgraph.file <- "./data/bed/ini-index.peak.pol-e-a.w1000.MA30-15.rep1.bedgraph"

output.name <- "ini-index.peak.pol-e-a.w1000.MA30-15.rep1.RT.bedgraph"
output.dir <- "."

# Replication-timing data derived from Zhao et al. (2020).
RT.wig <- "./data/wig/hct116.mean-RT.wig" 


# Import RT
RT.wig_file <- import(RT.wig)

RT.wig_file <- as.data.frame(RT.wig_file)
RT.wig_file <- RT.wig_file[, c(1, 2, 6)]
colnames(RT.wig_file) <- c("chro", "start", "RT")

# Import bedgraph
bg_file <- read.table(bedgraph.file, header = FALSE)
bg_file <- bg_file[, c(1,2, 4)]
colnames(bg_file) <- c("chro", "start", "score")

# Merge dataframes
RT.merged_df <- merge(RT.wig_file, bg_file, by = c("chro", "start"), all = TRUE)%>%
  arrange(chro, start)

# Function to fill missing RT values
fill_missing_rt <- function(df) {
  for (i in 1:nrow(df)) {
    if (is.na(df$RT[i])) {
      current_chro <- df$chro[i]
      start_pos <- df$start[i]
      
      # Find the closest non-NA RT value within the same chromosome
      non_na_indices <- which(df$chro == current_chro & !is.na(df$RT))
      
      # Check if there are any non-NA RT values within the same chromosome
      if (length(non_na_indices) > 0) {
        distances <- abs(df$start[non_na_indices] - start_pos)
        closest_index <- non_na_indices[which.min(distances)]
        
        # Fill the missing RT value
        df$RT[i] <- df$RT[closest_index]
      }
    }
  }
  return(df)
}

# Apply the function to fill missing RT values
RT.merged_df <- fill_missing_rt(RT.merged_df)

# Remove incomplete cases and process for export as bedgraph
bg_RT <- RT.merged_df[complete.cases(RT.merged_df), ]
bg_RT$end <- bg_RT$start + 1
bg_RT <- bg_RT[, c(1, 2, 5, 3)]

bedgraph.name <- file.path(output.dir, output.name)
if(!dir.exists(output.dir))dir.create(output.dir, recursive = TRUE)

write.table(
  bg_RT,
  file = bedgraph.name,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)
