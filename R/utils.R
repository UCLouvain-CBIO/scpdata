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
  if(nrow(fData(obj)) == 0) stop("'fData(obj)' cannot be empty")
  return(aggregateByProtein(obj))
}



scp_image <- function(obj, row.ord = NULL, col.ord = NULL){
  M <- exprs(obj)
  
  # Order the rows
  if(!is.null(row.ord)){
    if(!row.ord %in% colnames(fData(obj))) stop("'row.ord' should be a column name in 'fData(obj)'")
    M <- M[order(fData(obj)[, row.ord]), ]
  }
  
  # Order the columns 
  if(!is.null(col.ord)){
    if(!col.ord %in% colnames(pData(obj))) stop("'col.ord' should be a column name in 'pData(obj)'")
    M <- M[, order(pData(obj)[, col.ord])]
  }
  
  # Create the heatmap
  image(t(M), xlab = "Cell index", ylab = "Feature index", useRaster = TRUE,
        col = colorRampPalette(c("coral", "gold1", "#08306B"))(1000),
        axes = FALSE)
  
}

scp_plotStats <- function(obj, what = c("both", "cells", "features"), 
                          xstat = "mean", ystat = "sd"){
  # Check and format arguments 
  what <- match.arg(what, c("both", "cells", "features"))
  if(is.character(xstat)){
    xstat <- switch(xstat, 
                    mean = list(fun = function(x) mean(x, na.rm = TRUE), title = "Mean"),
                    median = list(fun = function(x) round(median(x, na.rm = TRUE), 15), title = "Median"),
                    sd = list(fun = function(x) sd(x, na.rm = TRUE), title = "Standard deviation"))
  }
  if(is.character(ystat)){
    ystat <- switch(ystat, 
                    mean = list(fun = function(x) mean(x, na.rm = TRUE), title = "Mean"),
                    median = list(fun = function(x) round(median(x, na.rm = TRUE), 15), title = "Median"),
                    sd = list(fun = function(x) sd(x, na.rm = TRUE), title = "Standard deviation"))
  }
  M <- exprs(obj)
  
  # Generate plots
  if(what %in% c("both", "features")){
    p1 <- ggplot(data = data.frame(var1 = apply(M, 1, xstat$fun),
                                   var2 = apply(M, 1, ystat$fun))) + 
      geom_point(aes(x = var1, y = var2), col = rgb(0, 0, 0.5, 0.5)) +
      xlab(xstat$title) + ylab(ystat$title) +
      ggtitle("Distribution statistics for features")
  }
  if(what %in% c("both", "cells")){
    p2 <- ggplot(data = data.frame(var1 = apply(M, 2, xstat$fun),
                                   var2 = apply(M, 2, ystat$fun))) + 
      geom_point(aes(x = var1, y = var2), col = rgb(0, 0, 0.5, 0.5)) +
      ggtitle("Distribution statistics for single cells") +
      xlab(xstat$title) + ylab(ystat$title)
  }
  
  # Print plots
  if(what == "both"){
    grid.arrange(p1, p2, nrow = 1)
  } else if (what == "features"){
    print(p1)
  } else if (what == "cells"){
    print(p2)
  }
}


scp_plotMissing <- function(obj, what = c("both", "cells", "features")){
  # Check and format arguments 
  what <- match.arg(what, c("both", "cells", "features")) 
  
  M <- exprs(obj)
  
  if(what %in% c("both", "features")){
    mis <- apply(M, 1, function(x) sum(is.na(x))/length(x))
    p1 <- ggplot(data = data.frame(mis), mapping = aes(x = mis)) + 
      geom_histogram(binwidth = 0.025, fill = "grey60", col = "grey40") +
      xlim(0, 1) + xlab("Missingness (%)") + ggtitle("Missing data among features")
  }
  if(what %in% c("both", "cells")){
    mis <- apply(M, 2, function(x) sum(is.na(x))/length(x))
    p2 <- ggplot(data = data.frame(mis), mapping = aes(x = mis)) + 
      geom_histogram(binwidth = 0.025, fill = "grey60", col = "grey40") +
      xlim(0, 1) + xlab("Missingness (%)") + ggtitle("Missing data among cells")
  }
  
  # Print plots
  if(what == "both"){
    grid.arrange(p1, p2, nrow = 1)
  } else if (what == "features"){
    print(p1)
  } else if (what == "cells"){
    print(p2)
  }
  
  
}


scp_plotCV <- function(obj, colorBy = NULL){
  dat <- exprs(obj)
  prots <- as.character(fData(obj)[,1])
  CVs <- do.call(rbind, lapply(unique(prots), function(prot){
    .idx <- prots == prot
    if(sum(.idx) <= 1) return(NULL)
    xx <- 2^dat[.idx,]
    CV <- apply(xx, 2, sd, na.rm = TRUE)/apply(xx, 2, mean, na.rm = TRUE)
    return(CV)
  }))
  dimnames(CVs) <- list(NULL, NULL)
  CVs <- melt(CVs)
  # CVs <- cbind(CVs, type = pData(obj)$celltype[CVs$Var2])
  ggplot(data = CVs, aes(x = Var2, y = value)) +
    geom_point(col = rgb(0.3,0.3,0.3,0.2)) + 
    stat_summary(aes(y = value,group=1), fun.y = median, colour = "red",
                 geom = "point",  size = 0.5, group = 1) +
    xlab("Cell index") + ylab("Coefficient of variation") +
    ggtitle("Protein CV distribution per cell (black = peptides; red = median)")
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

