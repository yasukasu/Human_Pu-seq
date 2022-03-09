
## v3 - eliminating bins without aligned reads in both strands

setwd("/home/yasukazu/R")
#library("seqinr")
source("./script/sub/YD.library.v2.R")
source("./script/sub/sub.Z.normalising.R")
source("./script/sub/sub.output.WigFile.Vstep.v4.R")

library("Rcpp")
sourceCpp("./script/rcpp/moving_ave.cpp")

RFD.wig.files.v2 <-function(path.Ld.f,  path.Ld.r, path.Lg.f, path.Lg.r, prefixs, file.g.size,  location.out, names.chro, w.size, MA, bin.count.T, sd.range){

  cat("\nprocessing data of", prefixs, "...\n")

  ## data from the tatal genme ##

  chro.sizes = read.table(file.g.size, row.names =1)[names.chro,1]
  genome.size = sum(as.numeric(chro.sizes))
  bin.num  = genome.size/w.size


  cat("Genome Size:", genome.size/10^6, "Mb (",genome.size, "bp)\n")

  ## inputting data from files ##

  Ld.f.df <- read.csv(path.Ld.f)
  Ld.r.df <- read.csv(path.Ld.r)
  Lg.f.df <- read.csv(path.Lg.f)
  Lg.r.df <- read.csv(path.Lg.r)

  #bin.size = Ld.f.df[,c("pos")][2]-Ld.f.df[,c("pos")][1]

  tatal.Ld.f = sum(Ld.f.df[,c("count")])
  tatal.Ld.r = sum(Ld.r.df[,c("count")])
  tatal.Lg.f = sum(Lg.f.df[,c("count")])
  tatal.Lg.r = sum(Lg.r.df[,c("count")])


  ave.bin.Ld.f = tatal.Ld.f/bin.num
  ave.bin.Ld.r = tatal.Ld.r/bin.num
  ave.bin.Lg.f = tatal.Lg.f/bin.num
  ave.bin.Lg.r = tatal.Lg.r/bin.num

  bin.size = w.size

  score.rfd.list <- list()
  position.list <- list()

  score.ld.rfd.list <- list()
  position.ld.list <- list()
  
  score.lg.rfd.list <- list()
  position.lg.list <- list()
  
  for(chromo in names.chro){

    cat(chromo, "...")

    Ld.f.df.chr = Ld.f.df[Ld.f.df$chro==chromo,]
    Ld.r.df.chr = Ld.r.df[Ld.r.df$chro==chromo,]
    Lg.f.df.chr = Lg.f.df[Lg.f.df$chro==chromo,]
    Lg.r.df.chr = Lg.r.df[Lg.r.df$chro==chromo,]

    if(nrow(Ld.f.df.chr)!=nrow(Ld.r.df.chr) &&
       nrow(Lg.f.df.chr)!=nrow(Lg.r.df.chr) &&
       nrow(Lg.f.df.chr)!=nrow(Ld.f.df.chr)){stop("The number of bins between f and r was not identical.\n")}
    

    Ld.f.count = Ld.f.df.chr[,c("count")]
    Ld.r.count = Ld.r.df.chr[,c("count")]
    Lg.f.count = Lg.f.df.chr[,c("count")]
    Lg.r.count = Lg.r.df.chr[,c("count")]

    pos = Ld.f.df.chr[,c("pos")]

    ### eliminating empty bins

    rmv.ld.vecs = (Ld.f.count<bin.count.T & Ld.r.count<bin.count.T)
    rmv.lg.vecs = (Lg.f.count<bin.count.T & Lg.r.count<bin.count.T)
    
    Ld.f.count[rmv.ld.vecs] = NA
    Ld.r.count[rmv.ld.vecs] = NA
    Lg.f.count[rmv.lg.vecs] = NA
    Lg.r.count[rmv.lg.vecs] = NA


    ### caluculation

    Ld.f.ech <- Ld.f.count/ave.bin.Ld.f
    Ld.r.ech <- Ld.r.count/ave.bin.Ld.r
    Lg.f.ech <- Lg.f.count/ave.bin.Lg.f
    Lg.r.ech <- Lg.r.count/ave.bin.Lg.r

    score.rfd     = (Ld.f.ech-Ld.r.ech-Lg.f.ech+Lg.r.ech)/(Ld.f.ech+Ld.r.ech+Lg.f.ech+Lg.r.ech)
    score.ld.rfd  = (Ld.f.ech-Ld.r.ech)/(Ld.f.ech+Ld.r.ech)
    score.lg.rfd  = (-Lg.f.ech+Lg.r.ech)/(Lg.f.ech+Lg.r.ech)
    
    ### moving average
    
    MA.1 = ifelse(length(pos)>MA, MA, length(pos))

    score.rfd.ma    <- rcpp_moving_ave(score.rfd, MA.1)
    score.ld.rfd.ma <- rcpp_moving_ave(score.ld.rfd, MA.1)
    score.lg.rfd.ma <- rcpp_moving_ave(score.lg.rfd, MA.1)
    
    rmv.vec.b = is.na(score.rfd.ma)
    
    score.rfd.ma.ex    <- round(score.rfd.ma[!rmv.vec.b], 4)
    score.ld.rfd.ma.ex <- round(score.ld.rfd.ma[!rmv.ld.vecs], 4)
    score.lg.rfd.ma.ex <- round(score.lg.rfd.ma[!rmv.lg.vecs], 4)
    
    pos.ex    = pos[!rmv.vec.b]
    pos.ld.ex = pos[!rmv.ld.vecs]
    pos.lg.ex = pos[!rmv.lg.vecs]
    
    score.rfd.list    <- c(score.rfd.list, list(score.rfd.ma.ex))
    score.ld.rfd.list <- c(score.ld.rfd.list, list(score.ld.rfd.ma.ex))
    score.lg.rfd.list <- c(score.lg.rfd.list, list(score.lg.rfd.ma.ex))
    
    position.list     <- c(position.list, list(pos.ex))
    position.ld.list  <- c(position.ld.list, list(pos.ld.ex))
    position.lg.list  <- c(position.lg.list, list(pos.lg.ex))
    
    cat(" done.\n")
  }

  names(score.rfd.list) = names.chro
  names(score.ld.rfd.list) = names.chro
  names(score.lg.rfd.list) = names.chro
  
  names(position.list) = names.chro
  names(position.ld.list) = names.chro
  names(position.lg.list) = names.chro

  path.wig = paste(location.out, "/", prefixs[1], "_RFD_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".wig", sep="")
  path.ld.wig = paste(location.out, "/", prefixs[2], "_RFD_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".wig", sep="")
  path.lg.wig = paste(location.out, "/", prefixs[3], "_RFD_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".wig", sep="")

  path.bw = paste(location.out, "/", prefixs[1], "_RFD_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".bw", sep="")
  path.ld.bw = paste(location.out, "/", prefixs[2], "_RFD_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".bw", sep="")
  path.lg.bw = paste(location.out, "/", prefixs[3], "_RFD_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".bw", sep="")
  
  name.rfd = paste(prefixs[1], "_RFD_w", bin.size, "_MA", MA, sep="")
  name.ld.rfd = paste(prefixs[2], "_RFD_w", bin.size, "_MA", MA, sep="")
  name.lg.rfd = paste(prefixs[3], "_RFD_w", bin.size, "_MA", MA, sep="")
  h.line = 0

  output.WigFile.Vstep.v4(name.rfd, path.wig, score.rfd.list, position.list, "200,28,89", h.line, y.rng=c(-1,1))
  output.WigFile.Vstep.v4(name.ld.rfd, path.ld.wig, score.ld.rfd.list, position.ld.list, "200,28,89", h.line, y.rng=c(-1,1))
  output.WigFile.Vstep.v4(name.lg.rfd, path.lg.wig, score.lg.rfd.list, position.lg.list, "200,28,89", h.line, y.rng=c(-1,1))
  
  system(paste("wigToBigWig", path.wig, file.g.size, path.bw, "-clip"))
  system(paste("wigToBigWig", path.ld.wig, file.g.size, path.ld.bw, "-clip"))
  system(paste("wigToBigWig", path.lg.wig, file.g.size, path.lg.bw, "-clip"))
  
  ## normalised RFD using sd values

  score.rfd.z.list     = Z.normalising.list(score.rfd.list)
  score.ld.rfd.z.list  = Z.normalising.list(score.ld.rfd.list)
  score.lg.rfd.z.list  = Z.normalising.list(score.lg.rfd.list)
  
  score.rfd.z2.list    = lapply(score.rfd.z.list, "+", sd.range)
  score.ld.rfd.z2.list = lapply(score.ld.rfd.z.list, "+", sd.range)
  score.lg.rfd.z2.list = lapply(score.lg.rfd.z.list, "+", sd.range)
  
  score.rfd.norm.list    = lapply(score.rfd.z2.list, "/", sd.range*2)
  score.ld.rfd.norm.list = lapply(score.ld.rfd.z2.list, "/", sd.range*2)
  score.lg.rfd.norm.list = lapply(score.lg.rfd.z2.list, "/", sd.range*2)
  
  score.rfd.norm.list     = lapply(score.rfd.norm.list, round, 4)
  score.ld.rfd.norm.list  = lapply(score.ld.rfd.norm.list, round, 4)
  score.lg.rfd.norm.list  = lapply(score.lg.rfd.norm.list, round, 4)
  
  
  path.nm.wig = paste(location.out, "/", prefixs[1], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".wig", sep="")
  path.nm.ld.wig = paste(location.out, "/", prefixs[2], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".wig", sep="")
  path.nm.lg.wig = paste(location.out, "/", prefixs[3], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".wig", sep="")
  
  path.nm.bw = paste(location.out, "/", prefixs[1], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".bw", sep="")
  path.nm.ld.bw = paste(location.out, "/", prefixs[2], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".bw", sep="")
  path.nm.lg.bw = paste(location.out, "/", prefixs[3], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, "_BT", bin.count.T, ".bw", sep="")
  
  name.nm.rfd = paste(prefixs[1], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, sep="")
  name.nm.ld.rfd = paste(prefixs[2], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, sep="")
  name.nm.lg.rfd = paste(prefixs[3], "_RFD_nm-0to1-", sd.range, "sd_w", bin.size, "_MA", MA, sep="")
 
   h.line = 0.5
  
  output.WigFile.Vstep.v4(name.nm.rfd   , path.nm.wig,    score.rfd.norm.list,    position.list,    "200,28,89", h.line, y.rng=c(0,1))
  output.WigFile.Vstep.v4(name.nm.ld.rfd, path.nm.ld.wig, score.ld.rfd.norm.list, position.ld.list, "200,28,89", h.line, y.rng=c(0,1))
  output.WigFile.Vstep.v4(name.nm.lg.rfd, path.nm.lg.wig, score.lg.rfd.norm.list, position.lg.list, "200,28,89", h.line, y.rng=c(0,1))
  
  system(paste("wigToBigWig", path.nm.wig,    file.g.size, path.nm.bw, "-clip"))
  system(paste("wigToBigWig", path.nm.ld.wig, file.g.size, path.nm.ld.bw, "-clip"))
  system(paste("wigToBigWig", path.nm.lg.wig, file.g.size, path.nm.lg.bw, "-clip"))

}


## refrence genome: fasta
# cat("Extracting fasta data ...\n")
# file.genome = "./data/genome/GRCh38/GRCh38.fasta"
file.g.size = "./data/genome/GRCh38/GRCh38.chrom.sizes"

## genome.chr.list = read.fasta(file = file.genome)
## names.chro = names(genome.chr.list)
## genome.size = length(unlist(genome.chr.list))

names.chro = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
               "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
               "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
#genome.size = 3095693983  ## wrong?? 15/1/2019
w.size = 1000
MA     =　30
bin.count.T = 5  # minmum reads to output bin data


cat("Chromosomes:", names.chro, "\n")
cat("Window Size:", w.size, "bp\n")

location.main = "./"

location.wig  = file.path(location.main, "wig")
location.out  = file.path(location.wig,  "")

if(!dir.exists(location.wig))dir.create(location.wig)
if(!dir.exists(location.out))dir.create(location.out)

path.Ld.f = "./data/pol-e-usage.watson.w1000.MA30.rep1.wig"
path.Ld.r = "./data/pol-e-usage.crick.w1000.MA30.rep1.wig" 
path.Lg.f = "./data/pol-a-usage.watson.w1000.MA30.rep1.wig"
path.Lg.r = "./data/pol-a-usage.crick.w1000.MA30.rep1.wig"  


path.lead.f = "./data/pol-e-usage.watson.w1000.MA30.rep1.wig"
path.lead.r = "./data/pol-e-usage.crick.w1000.MA30.rep1.wig" 
path.lagg.f = "./data/pol-a-usage.watson.w1000.MA30.rep1.wig"
path.lagg.r = "./data/pol-a-usage.crick.w1000.MA30.rep1.wig"  

prefixs = c("hct116-2-pol", "hct116-pol-e", "hct116-pol-a")

RFD.wig.files.v2(path.Ld.f, path.Ld.r, path.Lg.f, path.Lg.r, prefixs, file.g.size, location.out, names.chro, w.size, MA, bin.count.T, 3)





