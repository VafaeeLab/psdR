#' complexity function
#'
#' Returns a complexity score of the data, based on Fischer's Discriminant Ratio.
#' Higher score indicates lower complexity
#'
#' @param data scRNA count data - genes Vs cells (dim : g x c)
#' @param classes cell annotation i.e. the class of each of the cells in data
#' @return complexity score : higher score indicates lower complexity
#' @export
complexity <- function(data, classes = NA){
  if (is_empty(classes)) {
    return (NA)
  }

  if (anyNA(classes)) {
    data <- data[, !is.na(classes)]
    classes <- classes[!is.na(classes)]
  }

  uniq_classes <- unique(classes)
  mean_df <- data.frame(row.names = rownames(data))
  var_df <- data.frame(row.names = rownames(data))

  for (i in c(1:length(uniq_classes))) {
    c <- uniq_classes[i]
    mean_df <- cbind(mean_df, data.frame(apply(data[, classes == c], 1, mean)))
    var_df <- cbind(var_df, data.frame(apply(data[, classes == c], 1, var)))
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
