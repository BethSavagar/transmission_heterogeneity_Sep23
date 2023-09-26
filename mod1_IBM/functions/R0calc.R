R0calc <- function(vT,vE){
  Table <- as.data.frame(table(vT))
  Table[,1] <- as.numeric(as.character(Table[,1]))
  Table[,3] <- rep(0,nrow(Table))
  
  for(i in 1:nrow(Table)){
    Table[i,3] <- sum(vE[vT==Table[i,1]]) # sum probability of exposure to a unit with i contacts
  }
  sum(Table[,3]) # probabilities should equal 1
  
  m <- matrix(nrow=nrow(Table),ncol=nrow(Table))
  
  for(i in 1:nrow(m)){
    for(j in 1:ncol(m)){
      
      m[i,j] <- Table[i,1]*Table[j,3] 
      
    }
  }
  
  R0est <- abs(eigen(m)$values)[1]
  return(R0est)
}