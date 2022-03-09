#setwd("~/R")
#source("./script/sub/YD.library.v2.R")
source("./script/sub/sub.import.Wig.into.list.R")
source("./script/sub/sub.output.WigFile.Vstep.v3.R")
source("./script/sub/sub.Z.normalising.R")


### version to nomolseid as Z(x)-Z(0)


library("Rcpp")
#sourceCpp("./script/rcpp/bincount.cpp")
sourceCpp("./script/rcpp/moving_ave.cpp")

### parmeter

MA = 30   # the parameter for moving ave (2N+1)

bin.size = 1000

prefix ="Pol-e-a"

version ="rep1"

### data from the tatal genme

file.g.size = "./data/genome/GRCh38.chrom.sizes"
chro.sizes = read.table(file.g.size, row.names =1, col.names=c("chro","size"))

###  outputting location

location.main    = "./"
location.wig.in  = file.path(location.main, "wig")
location.wig.out = file.path(location.main, "wig")


if(!dir.exists(location.wig.in) )stop("The directory of input files does not exist.:", location.in)
if(!dir.exists(location.wig.out))dir.create(location.wig.out)

### inpputted wig files

path.lead.f = "./data/pol-e-usage.watson.w1000.MA30.rep1.wig"
path.lead.r = "./data/pol-e-usage.crick.w1000.MA30.rep1.wig" 
path.lagg.f = "./data/pol-a-usage.watson.w1000.MA30.rep1.wig"
path.lagg.r = "./data/pol-a-usage.crick.w1000.MA30.rep1.wig"  

fnames=basename(c(path.lead.f, path.lead.r, path.lagg.f, path.lagg.r))
f.str = strsplit(fnames, "-|_")
ori.MA= gsub("MA", "", sapply(f.str, function(x)x[grep("MA", x)]))

if(length(unique(ori.MA))!=1){stop("MA values of inputed wig files are not identical.")}

###  inputing wig data into list 

list.wig.lead.f = import.Wig.into.list(path.lead.f)
list.wig.lead.r = import.Wig.into.list(path.lead.r)
list.wig.lagg.f = import.Wig.into.list(path.lagg.f)
list.wig.lagg.r = import.Wig.into.list(path.lagg.r)

names.chro = names(list.wig.lead.f)
if(!identical(names.chro, names(list.wig.lead.r)) ||
   !identical(names.chro, names(list.wig.lagg.f)) ||
   !identical(names.chro, names(list.wig.lagg.r))){stop("The numbers of chro. in the wig files (f & r) are not identical.\n")}

#### function for ini index 

calc.ini.index = function(diff.vec){
  
  if(anyNA(diff.vec) || length(diff.vec)!=4){
    return(NA)
  }
  
  val = diff.vec[1]-diff.vec[2]-diff.vec[3]+diff.vec[4]
  
  if(diff.vec[1]>0 && diff.vec[2]<0 && diff.vec[3]<0 && diff.vec[4]>0){
    return(val)
  } else if(diff.vec[1]<0 && diff.vec[2]>0 && diff.vec[3]>0 && diff.vec[4]<0){
    return(val)    
  } else {
    return(NA)    
  }
  
}

### main
ini.index.list <- list()

for(chromo in names.chro){
  
  cat(chromo, "...")
  
  ## chromosome condinate
  chro.len = chro.sizes[chromo,]
  chro.N = ceiling(chro.len/bin.size)
  pos.all =seq(bin.size/2, by=bin.size, length=chro.N)
  
  ## wig file data
  
  wig.lead.f.chr = list.wig.lead.f[[chromo]]
  wig.lead.r.chr = list.wig.lead.r[[chromo]]
  wig.lagg.f.chr = list.wig.lagg.f[[chromo]]
  wig.lagg.r.chr = list.wig.lagg.r[[chromo]]
  
  pos.lead.chr = as.integer(names(wig.lead.f.chr))
  pos.lagg.chr = as.integer(names(wig.lagg.f.chr))
  
 
  
  if(length(pos.lead.chr)!=length(wig.lead.r.chr) ||
     length(pos.lagg.chr)!=length(wig.lagg.r.chr)) {stop(chromo,": The numbers of bins in the wig files (f & r) are not identical.\n")}
  
  ## arraging data of all coodinates
  
  lead.f.chr = rep(NA, chro.N)
  lead.r.chr = rep(NA, chro.N)
  lagg.f.chr = rep(NA, chro.N)
  lagg.r.chr = rep(NA, chro.N)
  
  names(lead.f.chr) = pos.all
  names(lead.r.chr) = pos.all 
  names(lagg.f.chr) = pos.all
  names(lagg.r.chr) = pos.all  
  
  
  lead.f.chr[as.character(pos.lead.chr)] = wig.lead.f.chr
  lead.r.chr[as.character(pos.lead.chr)] = wig.lead.r.chr
  lagg.f.chr[as.character(pos.lagg.chr)] = wig.lagg.f.chr
  lagg.r.chr[as.character(pos.lagg.chr)] = wig.lagg.r.chr  
  
  ## calculation
  
  lead.diff.f.chr = diff(lead.f.chr)
  lead.diff.r.chr = diff(lead.r.chr)
  lagg.diff.f.chr = diff(lagg.f.chr)
  lagg.diff.r.chr = diff(lagg.r.chr)  
  
  lead.diff.f.chr.ma <- rcpp_moving_ave(lead.diff.f.chr, MA)
  lead.diff.r.chr.ma <- rcpp_moving_ave(lead.diff.r.chr, MA)
  lagg.diff.f.chr.ma <- rcpp_moving_ave(lagg.diff.f.chr, MA)
  lagg.diff.r.chr.ma <- rcpp_moving_ave(lagg.diff.r.chr, MA)
  
  
  ini.index.chr = apply(cbind(lead.diff.f.chr.ma, lead.diff.r.chr.ma, lagg.diff.f.chr.ma, lagg.diff.r.chr.ma), 1,  calc.ini.index)
  
  ini.index.chr <- round(ini.index.chr, 4)
  
  pos.diff = (pos.all+bin.size/2)[-length(pos.all)]
  
  names(ini.index.chr) = pos.diff
  
  ini.index.chr.ex <- round(ini.index.chr[!is.na(ini.index.chr)], 4)
  
  ini.index.list <- c(ini.index.list, list(ini.index.chr.ex))
  
  cat(" done.\n")
  
  
}

names(ini.index.list) = names.chro

# nomarlisation

cat("Normalising to Z(x)-Z(0)...")

ini.index.z.list = Z.normalising.list(ini.index.list)
ini.index.z_0 = Z.value.extract(0, unlist(ini.index.list))

ini.index.z2.list = lapply(ini.index.z.list, "-", ini.index.z_0)
  
cat("done.\n") 




path.wig.s =  paste(location.wig.out, "/", prefix, ".ini-index_norm_z.w", bin.size, ".MA", ori.MA[1], "-", MA,".", version, ".wig", sep="")

path.bw.s  =  paste(location.wig.out, "/", prefix, ".ini-index_norm_z.w", bin.size, ".MA" ,ori.MA[1], "-", MA, ".", version, ".bw", sep="")



name.diff.s =    paste(prefix, "_ini-index_w", bin.size, "_MA", ori.MA[1], "-", MA, sep="")

disc.s=paste0("Analysis by Y.Daigaku; Source:", paste(fnames, collapse=" "))

h.line = 0

output.WigFile.Vstep.v3(name.diff.s, path.wig.s, ini.index.z2.list, bin.size, "116,60,117", h.line, y.rng=c(-0.15,0.15), description=disc.s)

system(paste("wigToBigWig", path.wig.s, file.g.size, path.bw.s, "-clip"))


