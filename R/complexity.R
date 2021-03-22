#' complexity function
#'
#' Returns a complexity score of the data, based on Fischer's Discriminant Ratio
#' Higher score indicates lower complexity
#'
#' @param data scRNA data in SingleCellExperiment object,
#'              with counts containing count data (genes Vs cells)
#'              and class information present in one of the colData
#' @param class_colname name of the column in colData containing class information
#' @return complexity score : higher score indicates lower complexity
#' @export
complexity <- function(data, class_colname){
  x <- as.matrix(SingleCellExperiment::counts(data))
  classes <- SingleCellExperiment::colData(data)[class_colname]
  classes <- classes[colnames(x),]

  x <- x[, !is.na(classes)]
  classes <- classes[!is.na(classes)]

  uniq_classes <- unique(classes)
  mean_df <- data.frame(row.names = rownames(x))
  var_df <- data.frame(row.names = rownames(x))

  for (i in c(1:length(uniq_classes))) {
    c <- uniq_classes[i]
    mean_df <- cbind(mean_df, data.frame(apply(x[, classes == c], 1, mean)))
    var_df <- cbind(var_df, data.frame(apply(x[, classes == c], 1, var)))
  }
  F <- c()
  for ( i in c(1 : (length(uniq_classes)-1) ) ) {
    for ( j in c((i+1) : length(uniq_classes)) ) {
      f_num <- (mean_df[, i] - mean_df[, j]) ^ 2
      f_den <- var_df[, i] + var_df[, j]
      f_den[f_den == 0] = .Machine$double.eps
      f <- f_num / f_den
      f <- caTools::trapz(c(1:length(f)), f)
      F <- c(F, f)
    }
  }
  complexity_score = mean(F)
  return (complexity_score)
}
