
## v3 - eliminating bins without aligned reads in both strands

source("./script/sub/sub.output.WigFile.Vstep.v3.R")
library("Rcpp")
sourceCpp("./script/rcpp/moving_ave.cpp")

rNTP.ratio.wig.files <-function(Strain, Control, file.g.size, location.in, location.out, chro_list, w.size, MA, bin.count.T){

  cat("\nprocessing data of", Strain, "...\n")

  ## data from the tatal genome ##

  chro.sizes = read.table(file.g.size, row.names =1)[chro_list,1]
  genome.size = sum(as.numeric(chro.sizes))
  bin.num  = genome.size/w.size


  cat("Genome Size:", genome.size/10^6, "Mb (",genome.size, "bp)\n")

  ## reading data from files ##

  fnames = list.files(location.in)

  st.fv.f = grep(paste0(Strain,".*\\.f-w", w.size, "\\.count.csv$"), fnames)
  st.fv.r = grep(paste0(Strain,".*\\.r-w", w.size, "\\.count.csv$"), fnames)
  ct.fv.f = grep(paste0(Control,".*\\.f-w", w.size, "\\.count.csv$"), fnames)
  ct.fv.r = grep(paste0(Control,".*\\.r-w", w.size, "\\.count.csv$"), fnames)

  if(length(st.fv.f)==0 || length(st.fv.r)==0){
    stop("The sample file was not detected in ", location.in)
  } else if(length(st.fv.f)>1 || length(st.fv.r)>1){
    stop("Multiple sample files are detected. f: ", length(st.fv.f), "r:", length(st.fv.r))
  }

  if(length(ct.fv.f)==0 || length(ct.fv.r)==0){
    stop("The cotrol file was not detected in ", location.in)
  } else if(length(ct.fv.f)>1 || length(ct.fv.r)>1){
    stop("Multiple control files are detected. f: ", length(ct.fv.f), "r:", length(ct.fv.r))
  }

  cat("Detected sample file:", fnames[st.fv.f], "&", fnames[st.fv.r], "\n")
  cat("Detected control file:", fnames[ct.fv.f], "&", fnames[ct.fv.r], "\n")

  path.mut.f = paste(location.in, "/", fnames[st.fv.f], sep="")
  path.mut.r = paste(location.in, "/", fnames[st.fv.r], sep="")
  path.ctl.f = paste(location.in, "/", fnames[ct.fv.f], sep="")
  path.ctl.r = paste(location.in, "/", fnames[ct.fv.r], sep="")


  mut.f.df <- read.csv(path.mut.f)
  mut.r.df <- read.csv(path.mut.r)
  ctl.f.df <- read.csv(path.ctl.f)
  ctl.r.df <- read.csv(path.ctl.r)

  #bin.size = mut.f.df[,c("pos")][2]-mut.f.df[,c("pos")][1]

  tatal.mut.f = sum(mut.f.df[,c("count")])
  tatal.mut.r = sum(mut.r.df[,c("count")])
  tatal.ctl.f = sum(ctl.f.df[,c("count")])
  tatal.ctl.r = sum(ctl.r.df[,c("count")])


  ave.bin.mut.f = tatal.mut.f/bin.num
  ave.bin.mut.r = tatal.mut.r/bin.num
  ave.bin.ctl.f = tatal.ctl.f/bin.num
  ave.bin.ctl.r = tatal.ctl.r/bin.num

  bin.size = w.size

  ratio.f.list <- list()
  ratio.r.list <- list()


  for(chromo in chro_list){

    cat(chromo, "...")

    mut.f.data.chr = mut.f.df[mut.f.df$chro==chromo,]
    mut.r.data.chr = mut.r.df[mut.r.df$chro==chromo,]
    ctl.f.data.chr = ctl.f.df[ctl.f.df$chro==chromo,]
    ctl.r.data.chr = ctl.r.df[ctl.r.df$chro==chromo,]

    if(nrow(mut.f.data.chr)!=nrow(mut.r.data.chr) &&
       nrow(ctl.f.data.chr)!=nrow(ctl.r.data.chr) &&
       nrow(ctl.f.data.chr)!=nrow(mut.f.data.chr)){stop("The number of bins between f and r was not identical.\n")}

    mut.f.count = mut.f.data.chr[,c("count")]
    mut.r.count = mut.r.data.chr[,c("count")]
    ctl.f.count = ctl.f.data.chr[,c("count")]
    ctl.r.count = ctl.r.data.chr[,c("count")]

    pos = mut.f.data.chr[,c("pos")]

    ### eliminating empty bins

    mut.0.vecs = (mut.f.count<bin.count.T & mut.r.count<bin.count.T)
    ctl.0.vecs = (ctl.f.count<bin.count.T | ctl.r.count<bin.count.T)

    rmv.vecs = (mut.0.vecs | ctl.0.vecs)

    mut.f.count[rmv.vecs]=NA
    mut.r.count[rmv.vecs]=NA
    ctl.f.count[rmv.vecs]=NA
    ctl.r.count[rmv.vecs]=NA


    ### caluculation

    mut.f.ech <- mut.f.count/ave.bin.mut.f
    mut.r.ech <- mut.r.count/ave.bin.mut.r
    ctl.f.ech <- ctl.f.count/ave.bin.ctl.f
    ctl.r.ech <- ctl.r.count/ave.bin.ctl.r


    ratio.f <- (mut.f.ech/ctl.f.ech)
    ratio.r <- (mut.r.ech/ctl.r.ech)


    ### moving average

    MA.1 = ifelse(length(pos)>MA, MA, length(pos))

    ratio.f.ma <- rcpp_moving_ave(ratio.f, MA.1)
    ratio.r.ma <- rcpp_moving_ave(ratio.r, MA.1)

    ratio.f.ma.ex <- round(ratio.f.ma[!rmv.vecs], 4)
    ratio.r.ma.ex <- round(ratio.r.ma[!rmv.vecs], 4)

    pos.ex = pos[!rmv.vecs]

    names(ratio.f.ma.ex)=pos.ex
    names(ratio.r.ma.ex)=pos.ex

    ratio.f.list <- c(ratio.f.list, list(ratio.f.ma.ex))
    ratio.r.list <- c(ratio.r.list, list(ratio.r.ma.ex))

    cat(" done.\n")
  }

  names(ratio.f.list) = chro_list
  names(ratio.r.list) = chro_list



  path.wig.f = paste(location.out, "/", Strain, "_pol-usage_watson_w", bin.size, "_MA", MA, ".wig", sep="")
  path.wig.r = paste(location.out, "/", Strain, "_pol-usage_crick_w",  bin.size, "_MA", MA, ".wig", sep="")

  path.bw.f = paste(location.out, "/", Strain, "_pol-usage_watson_w", bin.size, "_MA", MA, ".bw", sep="")
  path.bw.r = paste(location.out, "/", Strain, "_pol-usage_crick_w",  bin.size, "_MA", MA, ".bw", sep="")

  name.f =    paste(Strain, "_watson-str_w", bin.size, "_MA", MA, sep="")
  name.r =    paste(Strain, "_crick-str_w",  bin.size, "_MA", MA, sep="")


  h.line = 1

  output.WigFile.Vstep.v3(name.f, path.wig.f, ratio.f.list, bin.size, "200,28,89", h.line, y.rng=c(0,2.5))
  output.WigFile.Vstep.v3(name.r, path.wig.r, ratio.r.list, bin.size, "200,28,89", h.line, y.rng=c(0,2.5))

  chrom.sizes=""

  system(paste("wigToBigWig", path.wig.f, file.g.size, path.bw.f, "-clip"))
  system(paste("wigToBigWig", path.wig.r, file.g.size, path.bw.r, "-clip"))


}


## genome

file.g.size = "./data/genome/GRCh38.chrom.sizes"

names.chro = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
               "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
               "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")

w.size = 1000
MA     =　30
bin.count.T = 5  # minmum reads to output bin data

cat("Chromosomes:", names.chro, "\n")
cat("Window Size:", w.size, "bp\n")

location.main = "./data"
location.in   = file.path(location.main, "count.csv")
location.out  =  "./"

if(!dir.exists(location.in))stop("The directory of input files does not exist.:", location.in)

if(!dir.exists(location.out))dir.create(location.out)


rNTP.ratio.wig.files("pol-a-y865f" , "pol-a-plus", file.g.size, location.in, location.out, names.chro, w.size, MA, bin.count.T)
rNTP.ratio.wig.files("pol-e-m630f" , "pol-e-plus", file.g.size, location.in, location.out, names.chro, w.size, MA, bin.count.T)
