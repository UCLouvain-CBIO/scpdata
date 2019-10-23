#' Normalization of single-cell proteomics data
#'
#' The function \code{scp_normalise} (and identically \code{scp_normalize}) 
#' performs normalization of cells or peptides/proteins or both of quantified
#' expression data (\code{MSnSet} objects). Different normalization methods are 
#' implemented (see Details). 
#'
#' @param obj 
#' An MSnSet object
#' @param what 
#' A character indicating whether columns (\code{"col"}), rows (\code{"row"}) or
#' both (\code{"both"}) should be normalized
#' @param method 
#' The type of normalization (see Details)
#'
#' @details 
#' 
#' \code{method == 1}: each row (if \code{what == "row"}) or column (if 
#' \code{what == "col"}) is subtracted by summary of the corresponding values. 
#' The summary function is the mean for rows and the median for columns. 
#'
#' @return
#' 
#' An \code{MSnSet} object similar to the input object \code{obj}, but where the 
#' expression values have been replaced by normalized expression values. 
#' 
#' @export
#'
#' @examples
#' data(specht2019)
#' scpNorm <- scp_normalize(specht2019, what = "both")
#' 
scp_normalize <- scp_normalise <- function(obj, what = "col", method = 1){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' should be an MSnSet object." )
  what <- match.arg(what, choices = c("col", "row", "both"))
  
  ## Method 1: implemented by Specht et al. 2019
  if(method == 1){
    # Normalize rows
    if(what %in% c("row", "both")) obj <- rowNormalize(obj)
    # Normalize columns
    if(what %in% c("col", "both")) obj <- colNormalize(obj)
  } else {
    stop("Method '", method, "' not implemented.")
  }
  return(obj)
}


#' Aggregation of peptide expression data into protein expression data
#'
#' The \code{scp_AggregateByProtein} function takes an \code{MSnSet} object and 
#' merges the peptides (rows) of the expression data into proteins. This merging
#' is performed by taking the median value for each cell (column) for a given
#' protein group. The peptide to protein mapping should be present in the 
#' feature data of the \code{\link{MSnSet}} object.
#'
#' @param obj 
#' An \code{MSnSet} object containing the peptide expression data.
#'
#' @return
#' A new \code{MSnSet} object containing the aggregated protein expression data.
#' Note that the feature data is adapted and will contain only protein 
#' annotations.
#' 
#' @export
#'
#' @examples
#' data(specht2019)
#' scProt <- scp_aggregateByProtein(specht2019)
#' 
scp_aggregateByProtein <- function(obj){
  if(nrow(fData) == 0) stop("'fData(obj)' cannot be empty")
  return(aggregateByProtein(obj))
}



####---- SPECHT ET AL. 2019 FUNCTIONS ----####

# This part of the script contains the code/algorithms used by Specht et al 
# (2019) for processing their SCP data

rowNormalize <- function(obj){
  exprs(obj)  <- t(apply(exprs(obj), 1, function(x) x - mean(x, na.rm = TRUE)))
  return(obj)
}

colNormalize <- function(obj){
  exprs(obj) <- apply(exprs(obj), 2, function(x) x - median(x, na.rm = TRUE))
  return(obj)
}

aggregateByProtein <- function(obj){
  prots <- fData(obj)[,1]
  x <- do.call(rbind, lapply(unique(prots), function(prot){
    xx <- exprs(obj)[prots == prot, , drop = F]
    apply(xx, 2, median, na.rm = TRUE)
  }))
  obj.new <- MSnSet(exprs = x, 
                    fData = data.frame(protein = unique(prots)),
                    pData = pData(obj),
                    experimentData = experimentData(obj))
  return(obj.new)
}


#' Impute missing values using K-nearest neighbours
#' 
#' Internal function 
#' @export
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

batchCorrect <- function(obj, batch, target){
  if(is.character(batch)){
    batch <- pData(obj)[, batch]
  } else if (!is.factor(batch)){
    stop("'batch' should be either a column name (character) in pData(obj) or a factor")
  }
  if(is.character(target)){
    target <- model.matrix(~ as.factor(pData(obj)[, target]))
  } else if (!is.matrix(target)){
    stop("'target' should be either a column name (character) in pData(obj) or a design matrix")
  }
  exprs(obj) <- ComBat(dat = exprs(obj), batch = batch, mod = target, par.prior = T)
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
  ggplot(data = CVs, aes(x = Var2, y = value)) +
    geom_point(col = rgb(0.3,0.3,0.3,0.2)) + 
    stat_summary(aes(y = value,group=1), fun.y = median, colour = "red",
                 geom = "point",  size = 0.5, group = 1) +
    # geom_hline(yintercept = 0.43,  col = "grey40") +
    xlab("Cell index") + ylab("Coefficient of variation") +
    ggtitle("Protein CV distribution per cell (black = peptides; red = median)")
  
}

show_heatmap <- function(obj, xlab = "Peptide index", ylab = "Cell index", ...){
  dat <- exprs(obj)
  image(x = 1:nrow(dat), y = 1:ncol(dat), z = dat, xlab = xlab, ylab = ylab, ...)
}

