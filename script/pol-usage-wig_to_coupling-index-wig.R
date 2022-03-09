source("./script/sub/sub.import.Wig.into.list.v2.R")
source("./script/sub/sub.output.WigFile.Vstep.v4.R")
source("./script/sub/sub.Z.normalising.R")


### version to nomolseid as Z(x)-Z(0)


library("Rcpp")
sourceCpp("./script/rcpp/moving_ave.cpp")

### parmeter

bin.size = 1000

prefix ="Pol-e-a"

z.norm = "no"

yl.wig = c(-0.1,0.1)

### data from the tatal genme

file.g.size = "./data/genome/GRCh38.chrom.sizes"
chro.sizes = read.table(file.g.size, row.names =1, col.names=c("chro","size"))

###  outputting location

location.main    = "./"
location.wig.out = "./"


if(!dir.exists(location.wig.out))dir.create(location.wig.out)

### inpputted wig files

path.lead.f = "./data/wig/pol-e-m630f_pol-usage_watson_w1000_MA30.bw"
path.lead.r = "./data/wig/pol-e-m630f_pol-usage_crick_w1000_MA30.bw" 
path.lagg.f = "./data/wig/pol-a-y865f_pol-usage_watson_w1000_MA30.bw"
path.lagg.r = "./data/wig/pol-a-y865f_pol-usage_crick_w1000_MA30.bw"  

fnames=basename(c(path.lead.f, path.lead.r, path.lagg.f, path.lagg.r))
f.str = strsplit(fnames, "-|_|\\.")
ori.MA = gsub("MA", "", sapply(f.str, function(x)x[grep("MA", x)]))

if(length(unique(ori.MA))!=1){stop("MA values of inputed wig files are not identical.")}

###  inputing wig data into list 

list.wig.lead.f = import.Wig.into.list.v2(path.lead.f)
list.wig.lead.r = import.Wig.into.list.v2(path.lead.r)
list.wig.lagg.f = import.Wig.into.list.v2(path.lagg.f)
list.wig.lagg.r = import.Wig.into.list.v2(path.lagg.r)

names.chro = names(list.wig.lead.f[["value"]])
if(!identical(names.chro, names(list.wig.lead.r[["value"]])) ||
   !identical(names.chro, names(list.wig.lagg.f[["value"]])) ||
   !identical(names.chro, names(list.wig.lagg.r[["value"]]))){stop("The numbers of chro. in the wig files (f & r) are not identical.\n")}

#### function for ini index 

calc.CI.f.index = function(value){
  
  if(anyNA(value[c(1,4)]) || length(value)!=4){
    
    return(NA)
    
 # } else if(value[4]>1 || value[1]>1){
  } else {
    
    val = (value[1]-value[4])/(value[1]+value[4])
    return(val)
  
  # } else {
  #   return(NA)
  }

}



calc.CI.r.index = function(value){
  
  if(anyNA(value[c(2,3)]) || length(value)!=4){
    
    return(NA)
    
#  } else if(value[3]>1 || value[2]>1){
  } else {
  
    val = (value[2]-value[3])/(value[2]+value[3])
    return(val)
    
  # } else {
  #   return(NA)
  }
}


### main

ave.lead.f = mean(unlist(list.wig.lead.f[["value"]]))
ave.lead.r = mean(unlist(list.wig.lead.r[["value"]]))
ave.lagg.f = mean(unlist(list.wig.lagg.f[["value"]]))
ave.lagg.r = mean(unlist(list.wig.lagg.r[["value"]]))





pos.f.list <- c()
pos.r.list <- c()
pol.cpl.index.f.list <- list()
pol.cpl.index.r.list <- list()

for(chromo in names.chro){
  
  cat(chromo, "...")
  
  ## chromosome condinate
  chro.len = chro.sizes[chromo,]
  chro.N   = ceiling(chro.len/bin.size)
  pos.chr  = seq(bin.size/2, by=bin.size, length=chro.N)
  
  ## wig file data
  
  wig.lead.f.chr = list.wig.lead.f[["value"]][[chromo]]
  wig.lead.r.chr = list.wig.lead.r[["value"]][[chromo]]
  wig.lagg.f.chr = list.wig.lagg.f[["value"]][[chromo]]
  wig.lagg.r.chr = list.wig.lagg.r[["value"]][[chromo]]
  
  pos.lead.chr = list.wig.lead.f[["position"]][[chromo]]
  pos.lagg.chr = list.wig.lagg.f[["position"]][[chromo]]
  
  
  if(length(pos.lead.chr)!=length(wig.lead.r.chr) ||
     length(pos.lagg.chr)!=length(wig.lagg.r.chr)) {stop(chromo,": The numbers of bins in the wig files (f & r) are not identical.\n")}
  
  ## arraging data of all coodinates
  
  lead.f.chr = rep(NA, chro.N)
  lead.r.chr = rep(NA, chro.N)
  lagg.f.chr = rep(NA, chro.N)
  lagg.r.chr = rep(NA, chro.N)
  
  names(lead.f.chr) = pos.chr
  names(lead.r.chr) = pos.chr 
  names(lagg.f.chr) = pos.chr
  names(lagg.r.chr) = pos.chr  
  # 
  
  lead.f.chr[as.character(pos.lead.chr)] = wig.lead.f.chr/ave.lead.f
  lead.r.chr[as.character(pos.lead.chr)] = wig.lead.r.chr/ave.lead.r
  lagg.f.chr[as.character(pos.lagg.chr)] = wig.lagg.f.chr/ave.lagg.f
  lagg.r.chr[as.character(pos.lagg.chr)] = wig.lagg.r.chr/ave.lagg.r 
  
  ## calculation
  
  Pol.scores.mat = cbind(lead.f.chr, lead.r.chr, lagg.f.chr, lagg.r.chr)
  
  pol.cpl.index.f.chr = apply(Pol.scores.mat, 1,  calc.CI.f.index)
  pol.cpl.index.r.chr = apply(Pol.scores.mat, 1,  calc.CI.r.index)
  
  pol.cpl.index.f.list <- c(pol.cpl.index.f.list, list(pol.cpl.index.f.chr[!is.na(pol.cpl.index.f.chr)]))
  pol.cpl.index.r.list <- c(pol.cpl.index.r.list, list(pol.cpl.index.r.chr[!is.na(pol.cpl.index.r.chr)]))
  
  pos.f.list <- c(pos.f.list, list(pos.chr[!is.na(pol.cpl.index.f.chr)]))
  pos.r.list <- c(pos.r.list, list(pos.chr[!is.na(pol.cpl.index.r.chr)]))
  cat(" done.\n")

}

### normalisation using genome average

ave.all.f = mean(unlist(pol.cpl.index.f.list))
ave.all.r = mean(unlist(pol.cpl.index.r.list))

#pol.cpl.index.f.nm.list = lapply(pol.cpl.index.f.list, "/", ave.all.f)
#pol.cpl.index.r.nm.list = lapply(pol.cpl.index.r.list, "/", ave.all.r)


pol.cpl.index.f.nm.list <- lapply(pol.cpl.index.f.list, round, 4)
pol.cpl.index.r.nm.list <- lapply(pol.cpl.index.r.list, round, 4)



names(pol.cpl.index.f.nm.list) = names.chro
names(pol.cpl.index.r.nm.list) = names.chro

#### scaling; Z ####

if(z.norm=="yes" || z.norm=="y"){
  
  cat("Normalising to Z(x)-Z(1)...")
  
  pol.cpl.index.f.z.list = Z.normalising.list(pol.cpl.index.f.nm.list)
  pol.cpl.index.r.z.list = Z.normalising.list(pol.cpl.index.r.nm.list)
  
  pol.cpl.index.f.z_1 = Z.value.extract(1, unlist(pol.cpl.index.f.list))
  pol.cpl.index.r.z_1 = Z.value.extract(1, unlist(pol.cpl.index.r.list))
  
  pol.cpl.index.f.z2.list = lapply(pol.cpl.index.f.z.list, "-", pol.cpl.index.f.z_1)
  pol.cpl.index.r.z2.list = lapply(pol.cpl.index.r.z.list, "-", pol.cpl.index.r.z_1)
  
  pol.cpl.index.f.out.list = lapply(pol.cpl.index.f.z2.list, round, 4)
  pol.cpl.index.r.out.list = lapply(pol.cpl.index.r.z2.list, round, 4)
  
  cat("done.\n") 
  
  path.f.wig =  paste(location.wig.out, "/", prefix, ".coupling-index-f_norm_z.w", bin.size, ".MA", ori.MA[1], ".wig", sep="")
  path.r.wig =  paste(location.wig.out, "/", prefix, ".coupling-index-r_norm_z.w", bin.size, ".MA", ori.MA[1], ".wig", sep="")
  
  path.f.bw  =  paste(location.wig.out, "/", prefix, ".coupling-index-f_norm_z.w", bin.size, ".MA" ,ori.MA[1], ".bw", sep="")
  path.r.bw  =  paste(location.wig.out, "/", prefix, ".coupling-index-r_norm_z.w", bin.size, ".MA" ,ori.MA[1], ".bw", sep="")
  
} else {
  
  pol.cpl.index.f.out.list = pol.cpl.index.f.nm.list
  pol.cpl.index.r.out.list = pol.cpl.index.r.nm.list 
  
  path.f.wig =  paste(location.wig.out, "/", prefix, ".coupling-index.rightward.w", bin.size, ".MA", ori.MA[1], ".wig", sep="")
  path.r.wig =  paste(location.wig.out, "/", prefix, ".coupling-index.leftward.w", bin.size, ".MA", ori.MA[1], ".wig", sep="")
  
  path.f.bw  =  paste(location.wig.out, "/", prefix, ".coupling-index.rightward.w", bin.size, ".MA" ,ori.MA[1], ".bw", sep="")
  path.r.bw  =  paste(location.wig.out, "/", prefix, ".coupling-index.leftward.w", bin.size, ".MA" ,ori.MA[1], ".bw", sep="")
  
  
}

name.f = paste(prefix, "_coupling-index-f_w", bin.size, "_MA", ori.MA[1], sep="")
name.r = paste(prefix, "_coupling-index-r_w", bin.size, "_MA", ori.MA[1], sep="")



disc.s=paste0("Analysis by Y.Daigaku; Source:", paste(fnames, collapse=" "))

h.line = 0

output.WigFile.Vstep.v4(name.f, path.f.wig, pol.cpl.index.f.out.list, pos.f.list, bin.size, "116,60,117", h.line, y.rng=yl.wig, description=disc.s)
output.WigFile.Vstep.v4(name.r, path.r.wig, pol.cpl.index.r.out.list, pos.r.list, bin.size, "116,60,117", h.line, y.rng=yl.wig, description=disc.s)

system(paste("wigToBigWig", path.f.wig, file.g.size, path.f.bw, "-clip"))
system(paste("wigToBigWig", path.r.wig, file.g.size, path.r.bw, "-clip"))


