# Data in wig or Bigwig files will be inputted in list. The position of data values stores as their 'names'
# rtracklayer is needed to be loaded; library("rtracklayer")
#
# This function is compatible for both variable step and fixted step

suppressMessages(library("rtracklayer"))

import.Wig.into.list <- function(path){
  
  if(!file.exists(path)){stop("  import.wig.into list: '", path, "' does not exsit.\n")}
  
  bwf <- WIGFile(path)
  
  if(regexpr('\\.wig$', path, ignore.case = T)>0){
    
    cat("  import.wig.into list:", "WigFile >", path, "\n")
    bwf <- WIGFile(path)
    
  } else if(regexpr('\\.bigwig$|\\.bw$', path, ignore.case = T)>0){
    
    cat("  import.wig.into list:", "BigWigFile >", path, "\n")
    bwf <- BigWigFile(path)
    
  } else {
    stop("  import.wig.into list:", "File type error!")
  }
  
  track <-import(bwf)
  
  chros <-    track @ seqinfo @ seqnames 
  lengths <-  track @ seqnames @ lengths
  pos   <- track @ ranges @ start
  score <- track@ elementMetadata @listData $ score
  
  val.list <- list()
  for(i in 1:length(chros)){
    #for(i in 1:5){
    
    N=lengths[i]
    chro=chros[i]
    
    #cat("list.BigWig.data >> processing ", chro, " (N=", N, ")...\n", sep="")
    v = score[1:N]
    names(v) = pos[1:N]
    
    val.list = c(val.list, list(v))
    
    score = score[-(1:N)]
    pos = pos[-(1:N)]
    
  }
  
  names(val.list)=chros
  (val.list)
  
}



# Data in 'variablestep' wig files will be inputted in list containing multiple matrixs(row1:position, row2:score)
# rtracklayer is needed to be loaded; library("rtracklayer")

import.V.Wig.into.list <- function(path){
  
  if(!file.exists(path)){stop("  import.wig.into list: '", path, "' does not exsit.\n")}
  
  bwf <- WIGFile(path)
  track <-import(bwf)

  val.list <- list()

  for(i in 1:length(track @ listData)){


    chro  <- track @ listData [i] $`R Track` @ seqinfo @ seqnames
    start <- track @ listData [i] $`R Track` @ ranges @ start
    width <- track @ listData [i] $`R Track` @ ranges @ width

    pos = start+as.integer(width/2)
    score = track @ listData [i] $`R Track` $score

    if(length(pos)!=length(score))stop("lengths of pos and score are different.")


    data.mat = rbind(pos, score)
    rownames(data.mat) = c("position", "score")

    #names(score) = pos

    val.list[i] = list(data.mat)
    names(val.list)[i] = chro
  }

  (val.list)


}

normalise.list.data <- function(list.data){
  
  all.data = unlist(list.data)
  
  f.norm = function(x)(x-mean(all.data))/sd(all.data)
  
  list.norm.data = lapply(list.data, f.norm)
  
  
  return(list.norm.data)

}

