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
  data <- data[rowSums(data) != 0, ]  #filter out genes with 0 count
  methods <- append('none', methods)
  complexity_df <- data.frame(Method = character(), PSD = logical(),
                              Complexity = double(), TimeTaken = double())

  for (method in methods) {

    if (method == 'none') {
      pp_data <- data
    }
    else if (method == 'CPM') {
      pp_data <- edgeR::cpm(data)
    }
    else if (method == 'Linnorm') {
      pp_data <- Linnorm::Linnorm(data)
    }
    else {
      print(paste('Skipping unknown method', method))
      next
    }

    complexity_df <- update_complexity_df(complexity_df, pp_data, classes, method, FALSE)
    complexity_df <- update_complexity_df(complexity_df, pp_data, classes, method, TRUE)

    print('finished')
    print(method)
    print('*****************')
  }

  return (complexity_df)
}


update_complexity_df <- function(complexity_df, pp_data, classes, method, is_psd){
  start <- Sys.time()
  if (is_psd) {
    pp_data <- psdR::psd(pp_data, filter = FALSE)
                  #genes with 0 count are already filtered out
                  #       in the body of 'compareMethods' function
                  #       So passing filter=FALSE here
  }
  complexity <- psdR::complexity(pp_data, classes)
  end <- Sys.time()
  complexity_df <- rbind(complexity_df,
                         data.frame(Method = method, PSD = is_psd, Complexity = complexity,
                                    TimeTaken = difftime(end, start, units = 'secs')))
  return (complexity_df)
}


