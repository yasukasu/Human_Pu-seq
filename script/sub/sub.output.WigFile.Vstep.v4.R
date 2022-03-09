#output.WigFile.v4(name <Name appearing in the browser>, out.path <location of outputfiles>,
#               data.chr.list <list of data in each chomosome, names of each element => position, names of list => chromosome>
#               bin.size, COlor<like "0,0,255", h.line <pos. of vertical line>)
# This function is usable to any organism.

# v4 avoiding names() to retain postion data to speed up process

output.WigFile.Vstep.v4 <- function(name, out.path, data.chr.list,  pos.chr.list, color, h.line=0, description="Analysis by Y. Daigaku", y.rng=c(0,5), auto.scale="off"){
 
  options(scipen=6)
  
   
  cat("  output.WigFile.Vstep.v4: outputting wig file -", out.path, "\n")
  header1 = paste('track type=wiggle_0 name="', name,
                  '" description="', description, ';', date(), 
                  '" visibility=full autoScale=', auto.scale, 
                  ' color=', color,
                  ' yLineOnOff=on yLineMark=', h.line, 
                  ' viewLimits=', y.rng[1], ':', y.rng[2],
                  ' priority=10', sep="")
  
  N.chr <- length(data.chr.list)
  chros = names(data.chr.list)
  bin.size = min(diff(pos.chr.list[[1]]))
  
  if(!length(chros))stop("  output.WigFile.Vstep.v4: Names of chromosomes were not inputted as names of list elements.\n")
  
  cat("  output.WigFile.Vstep.v4: data of", N.chr, "chromosomes is being processed.\n")
  
  if(length(header1)!=1){cat("  output.WigFile.Vstep.v4: outputting wig file - Error. Data length of inpputted parameter is more than 2."); q(1);}
  
  out <- file(out.path, "w") 
  on.exit(close(out))
  
  writeLines(header1, out, sep="\n")
  
  for(i in 1:N.chr){
    
    score = data.chr.list[[i]]
    pos = pos.chr.list[[i]]
    
    if(length(score)!=length(pos))stop("  output.WigFile.Vstep.v4: ", names(data.chr.list)[i], " data length error!")
    
    Nc <- length(data.chr.list[[i]])
    
     # outputting data 
    # header2 = paste('variableStep chrom=', chros[i],
    #                 ' span=', bin.size, sep="")
    header2 = paste('variableStep chrom=', chros[i], sep="")
    writeLines(header2, out, sep="\n") 
    
    if(Nc==0)next
    
    for (ii in 1:Nc) {
      writeLines(paste(pos[ii], score[ii], sep=" "), out, sep="\n")
    }
    cat("  output.WigFile.Vstep.v4: ", chros[i], " => ", Nc, " data points (", min(pos), "-", max(pos),")\n", sep="")
    
  }
 
}