#' compare_methods function
#'
#' Returns a comparison of complexity scores of the data before and after
#' applying preprocessing methods and also after applying psdR.
#' Also plots the datasets in all above cases after applying tSNE
#'
#' @param data scRNA count data - genes Vs cells (dim : g x c)
#' @param classes cell annotation i.e. the class of each of the cells in the data
#' @param methods vector of preprocessing methods to be compared
#' @param psd logical value to indicate if psd should be applied while comparing methods
#' @param show_plots logical value to indicate if plots should be shown
#' @param perplexity perplexity value while computing tSNE
#' @param plot_file_name name of the file to save the plot : defaults to "tsne_plot.png"
#' @param plot_colour_label title of the plot legend
#' @param plot_title title of the plot
#' @return dataframe containing comparison of complexity scores
#' @export
compare_methods <- function(data, classes = NA, methods = c('CPM', 'Linnorm'),
                       psd = TRUE, show_plots = TRUE,
                       perplexity = 30, plot_file_name = "tsne_plot.png",
                       plot_colour_label = "Cell Types",
                       plot_title = "tSNE embeddings") {
  if (!is_empty(classes) && anyNA(classes)) {
    data <- data[, !is.na(classes)]
    classes <- classes[!is.na(classes)]
  }
  data <- data[rowSums(data) != 0, ]  #filter out genes with 0 count
  methods <- append('none', methods)
  df_list <- init_df_list(classes)
  allowed_methods <- c('none', 'CPM', 'NormCPM',
                       'TMM', 'Linnorm', 'Seurat')
  for (method in methods) {
    if (method %in% allowed_methods) {
      pp_data <- preprocess_data(data, method)
    }
    else {
      print(paste('Skipping unknown method', method))
      next
    }
    df_list <- update_df_list(df_list, pp_data, classes, method,
                                          FALSE, show_plots, perplexity)
    if (psd) {
      df_list <- update_df_list(df_list, pp_data, classes, method,
                                TRUE, show_plots, perplexity)
    }
    print(paste("Completed", method))
  }
  if (show_plots) {
    create_plot(df_list, methods, psd, plot_file_name, plot_colour_label,
                plot_title)
  }
  return (df_list[[1]])
}


preprocess_data <- function(data, method){
  if (method == 'none') {
    pp_data <- data
  }
  else if (method == 'CPM') {
    pp_data <- edgeR::cpm(data)
  }
  else if (method == 'NormCPM') {
    pp_data <- scale(edgeR::cpm(data))
  }
  else if (method == 'TMM') {
    pp_data <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(data), method = 'TMM'))
  }
  else if (method == 'Linnorm') {
    pp_data <- Linnorm::Linnorm(data)
  }
  else if (method == 'Seurat') {
    pp_data <- as.matrix(Seurat::NormalizeData(data))
  }
  else {
    pp_data <- data
  }
  return (pp_data)
}


init_df_list <- function(classes){
  complexity_df <- data.frame(Method = character(), Complexity = double(), TimeTaken = double())
  if (is_empty(classes)) {
    tsne_result_df <- data.frame(Method = character(), x = double(), y = double())
  }
  else {
    tsne_result_df <- data.frame(Method = character(), x = double(), y = double(),
                                 Colour = character())
  }
  return (list(complexity_df, tsne_result_df))
}


update_df_list <- function(df_list, pp_data, classes, method,
                                 is_psd, show_plots = TRUE,
                                perplexity, random_seed = 1){
  complexity_df <- df_list[[1]]
  tsne_result_df <- df_list[[2]]
  start <- Sys.time()
  if (is_psd) {
    pp_data <- psdR::psd(pp_data, filter = FALSE)
                  #genes with 0 count are already filtered out
                  #       in the body of 'compareMethods' function
                  #       So passing filter=FALSE here
    method <- paste(method, "with psd", sep = " - ")
  }
  else {
    method <- paste(method, "no psd", sep = " - ")
  }
  complexity <- psdR::complexity(pp_data, classes)
  end <- Sys.time()
  complexity_df <- rbind(complexity_df,
                         data.frame(Method = method, Complexity = complexity,
                                    TimeTaken = difftime(end, start, units = 'secs')))
  if (show_plots) {
    set.seed(random_seed)
    tsne_result <- Rtsne::Rtsne(t(pp_data), perplexity = perplexity)
    if (is_empty(classes)) {
      tsne_result_row <- data.frame(Method = method, x = tsne_result$Y[,1],
                                    y = tsne_result$Y[,2])
    }
    else {
      tsne_result_row <- data.frame(Method = method, x = tsne_result$Y[,1],
                                    y = tsne_result$Y[,2], Colour = classes)
    }
    tsne_result_df <- rbind(tsne_result_df, tsne_result_row)
  }
  return (list(complexity_df, tsne_result_df))
}


create_plot <- function(df_list, methods, psd, plot_file_name,
                        plot_colour_label,
                        plot_title){
  nrow <- length(methods)
  ncol <- if (psd) 2 else 1
  if ('Colour' %in% colnames(df_list[[2]])) {
    tsne_plot <- ggplot2::ggplot(df_list[[2]]) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = Colour)) +
      ggplot2::facet_wrap(ggplot2::vars(Method), nrow = nrow, ncol = ncol) +
      ggplot2::labs(colour = plot_colour_label, title = plot_title) +
      ggplot2::xlab("Dimension 1") +
      ggplot2::ylab("Dimension 2")
  }
  else{
    tsne_plot <- ggplot2::ggplot(df_list[[2]]) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y)) +
      ggplot2::facet_wrap(ggplot2::vars(Method), nrow = nrow, ncol = ncol) +
      ggplot2::labs(title = plot_title) +
      ggplot2::xlab("Dimension 1") +
      ggplot2::ylab("Dimension 2")
  }
  ggplot2::ggsave(plot_file_name, tsne_plot)
}


is_empty <- function(info) {
  return (length(info) == 0 || (length(info) == 1 && is.na(info)))
}
