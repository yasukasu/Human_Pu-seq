#output.WigFile.v3(name <Name appearing in the browser>, out.path <location of outputfiles>,
#               data.chr.list <list of data in each chomosome, names of each element => position, names of list => chromosome>
#               bin.size, COlor<like "0,0,255", h.line <pos. of vertical line>)
# This function is usable to any organism.



output.WigFile.Vstep.v3 <- function(name, out.path, data.chr.list, bin.size, color, h.line, description="Analysis by Y. Daigaku", y.rng=c(0,5), auto.scale="off"){
 
 
   
  cat("### Wig file -", out.path, "###\n")
  header1 = paste('track type=wiggle_0 name="', name,
                  '" description="', description, ';', date(), 
                  '" visibility=full autoScale=', auto.scale, 
                  ' color=', color,
                  ' yLineOnOff=on yLineMark=', h.line, 
                  ' viewLimits=', y.rng[1], ':', y.rng[2],
                  ' priority=10', sep="")
  
  N.chr <- length(data.chr.list)
  chros = names(data.chr.list)
  
  if(!length(chros))stop("Names of chromosomes were not inputted as names of list elements.\n")
  
  cat("Data of", N.chr, "chromosomes is being processed.\n")
  
  out <- file(out.path, "w") 
  on.exit(close(out))
  
  writeLines(header1, out, sep="\n")
  
  for(i in 1:N.chr){
    
    score = data.chr.list[[i]]
    pos = as.integer(names(score))
    
    Nc <- length(data.chr.list[[i]])
    
     # outputting data 
    header2 = paste('variableStep chrom=', chros[i],
                    ' span=', bin.size, sep="")
    writeLines(header2, out, sep="\n") 
    
    if(Nc==0)next
    
    for (ii in 1:Nc) {
      writeLines(paste(pos[ii], score[ii], sep=" "), out, sep="\n")
    }
    cat(chros[i], " => ", Nc, " bins (0-", bin.size*Nc,")\n", sep="")
    
  }
 
}