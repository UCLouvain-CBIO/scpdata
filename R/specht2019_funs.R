####################################################
#### SPECHT ET AL. 2019 DATA SPECIFIC FUNCTIONS ####
####################################################

row.normalize <- function(obj){
  exprs(obj)  <- t(apply(exprs(obj), 1, function(x) x - mean(x, na.rm = TRUE)))
  return(obj)
}

col.normalize <- function(obj){
  exprs(obj) <- apply(exprs(obj), 2, function(x) x - median(x, na.rm = TRUE))
  return(obj)
}

aggregateByProtein <- function(obj){
  prots <- fData(obj)[,1]
  x <- do.call(rbind, lapply(unique(prots), function(prot){
    xx <- exprs(obj)[prots == prot, , drop = F]
    apply(xx, 2, median, na.rm = TRUE)
  }))
  return(MSnSet(exprs = x, 
                fData = data.frame(protein = unique(prots)),
                experimentData = experimentData(obj)))
}

imputeKNN <- function(obj, k = 3){
  dat <- exprs(obj)
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
      }
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ])
      }
    }
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
  }
  exprs(obj) <- dat.imp
  return(obj)
}

plotCV <- function(obj){
  dat <- exprs(obj)
  prots <- as.character(fData(obj)[,1])
  CVs <- do.call(rbind, lapply(unique(prots), function(prot){
    .idx <- prots == prot
    if(sum(.idx) <= 1) return(NULL)
    xx <- 2^dat[.idx,]
    CV <- apply(xx, 2, sd, na.rm = TRUE)/apply(xx, 2, mean, na.rm = TRUE)
    return(CV)
  }))
  dimnames(CVs) <- list(NULL, colnames(dat))
  meds <- apply(CVs, 1, median, na.rm = TRUE)
  CVs <- melt(CVs)
  CVs$Var2 <- as.numeric(gsub(CVs$Var2, pattern = "[a-z]", replacement = ""))
  ggplot(data = CVs) +
    geom_point(aes(x = Var2, y = value), col = rgb(0.3,0.3,0.3,0.2)) + 
    geom_line(aes(x = unique(Var2), y = meds), col = "firebrick") +
    xlab("Cell index") + ylab("Coefficient of variation")
  
}

