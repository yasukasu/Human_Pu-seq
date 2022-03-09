# Data in wig or Bigwig files will be inputted in list. The position of data values stores as their 'names'
# rtracklayer is needed to be loaded; library("rtracklayer")
#
# This function is compatible for both variable step and fixted step
# v2 avoiding names() to retain postion data to speed up process

suppressMessages(library("rtracklayer"))

import.Wig.into.list.v2 <- function(path){
  
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
  pos.list <- list()
  
  for(i in 1:length(chros)){
    #for(i in 1:5){
    
    N=lengths[i]
    chro=chros[i]
    
    #cat("list.BigWig.data >> processing ", chro, " (N=", N, ")...\n", sep="")
    
  
    val.list = c(val.list, list(score[1:N]))
    pos.list = c(pos.list, list(pos[1:N]))
    
    score = score[-(1:N)]
    pos = pos[-(1:N)]
    
  }
  
  names(val.list)=chros
  names(pos.list)=chros
  
  list(value=val.list, position=pos.list)
  
  
  
}

