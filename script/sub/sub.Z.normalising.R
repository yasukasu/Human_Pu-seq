

Z.normalising.list <- function(list.data){
  
  list.N = length(list.data)
  
  data.all = unlist(list.data)
  data.len = sapply(list.data, length)
  
  ave = mean(data.all[!is.na(data.all)])
  sd  = sd(data.all[!is.na(data.all)])
  
  data.z.all = (data.all-ave)/sd
  
  list.Z = list()
  
  
  for(i in 1:list.N){
    
    if(!data.len[i]){
      list.Z[[i]]=integer(0)
      next
    }
    
    vv = 1:data.len[i]
    
    list.Z[[i]] = data.z.all[vv]
    
    data.z.all = data.z.all[-(vv)]
    
    names(list.Z[[i]]) = names(list.data[[i]])
    
  }
  
  names(list.Z) = names(list.data)
  
  (list.Z)
}


Z.value.extract <- function(x, vec.dist){

  ave = mean(vec.dist[!is.na(vec.dist)])
  sd  = sd(vec.dist[!is.na(vec.dist)])
  
  z=(x-ave)/sd
  
  (z)
}