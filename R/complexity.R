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
  # data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
  # class_colname <- "cell_type1"
  x <- as.matrix(counts(data))
  classes <- colData(data)[class_colname]
  classes <- classes[colnames(x),]
  classes <- unique(classes)

  complexity_score <- 0 #logic to be added
  return (complexity_score)
}
