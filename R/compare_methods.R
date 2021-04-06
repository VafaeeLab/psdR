#' compare_methods function
#'
#' Returns a comparison of complexity scores of the data before and after
#' applying preprocessing methods and also after applying psdR.
#' Also plots the datasets in all above cases after applying tSNE
#'
#' @param data scRNA count data - genes Vs cells (dim : g x c)
#' @param classes cell annotation i.e. the class of each of the cells in data
#' @param methods list of preprocessing methods to be compared
#' @param psd should psd be applied while comparing methods
#' @param show_plots  should plots be shown
#' @return comparison of complexity scores
#' @export
compare_methods <- function(data, classes = NA, methods = c('CPM', 'Linnorm'),
                       psd = TRUE, show_plots = TRUE) {
  if (anyNA(classes)) {
    data <- data[, !is.na(classes)]
    classes <- classes[!is.na(classes)]
  }
  data <- data[rowSums(data) > 0, ]  #filter out genes with 0 count
  methods <- append('none', methods)
  complexity_df <- data.frame(Method = character(), PSD = logical(),
                              Complexity = double(), TimeTaken = double())
  for (method in methods) {
    if (method == 'none') {
      pp_data <- data
    }
    else {
      if (method == 'CPM') {
        pp_data <- edgeR::cpm(data)
      }
      else if (method == 'Linnorm') {
        pp_data <- Linnorm::Linnorm(data)
      }
      else {
        print('Unknown method - no preprocesseing done')
        pp_data <- data
      }
    }

    start <- Sys.time()
    complexity <- psdR::complexity(pp_data, classes)
    end <- Sys.time()
    complexity_df <- rbind(complexity_df,
                           data.frame(Method = method, PSD = FALSE, Complexity = complexity,
                                      TimeTaken = difftime(end, start, units = 'secs')))

    start <- Sys.time()
    complexity <- psdR::complexity(psdR::psd(pp_data), classes)
    end <- Sys.time()
    complexity_df <- rbind(complexity_df,
                           data.frame(Method = method, PSD = TRUE, Complexity = complexity,
                                      TimeTaken = difftime(end, start, units = 'secs')))
    print('finished')
    print(method)
    print('*****************')
  }

  return (complexity_df)
}


data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
class_colname <- "cell_type1"
x <- as.matrix(SingleCellExperiment::counts(data))
classes <- SingleCellExperiment::colData(data)[class_colname]
classes <- classes[colnames(x),]

complexity_df <- compare_methods(x, classes)

# data <- x[1:3000, 1:300]
# classes <- classes[1:300]
# data <- x
#
# start <- Sys.time()
# if (anyNA(classes)) {
#   data <- data[, !is.na(classes)]
#   classes <- classes[!is.na(classes)]
# }
# data <- data[rowSums(data) > 0, ]  #filter out genes with 0 count
#
# method <- 'Linnorm'
# if (method == 'none') {
#   pp_data <- data
# } else {
#   if (method == 'CPM') {
#     pp_data <- edgeR::cpm(data)
#   }
#   else if (method == 'Linnorm') {
#     pp_data <- Linnorm::Linnorm(data)
#   }
#   else {
#     pp_data <- data
#   }
# }
#
# complexity <- psdR::complexity(pp_data, classes)
# complexity_df <- data.frame(Method = method, PSD = FALSE, Complexity = complexity)
# complexity_df <- rbind(complexity_df,
#                        data.frame(Method = method, PSD = FALSE, Complexity = complexity))
#
# complexity <- psdR::complexity(psdR::psd(pp_data), classes)
# complexity_df <- rbind(complexity_df,
#                        data.frame(Method = method, PSD = TRUE, Complexity = complexity))

start <- Sys.time()
end <- Sys.time()
data.frame(TimeTaken = difftime(end, start, units = 'secs'))

print(paste('execution time :', difftime(end, start, units = 'secs')))



