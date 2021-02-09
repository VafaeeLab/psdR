#' psd function
#'
#' Compute the Power Spectral Density (psd) transformation of the given scRNA count data
#'
#' @param data scRNA count data - cells Vs genes
#' @return psd transformation of input data
#' @export
psd <- function(data){
  data.cor <- cor(data)         #dim : g x g
  data.psd <- data %*% data.cor #dim : c x g
  #do remaining steps of psd
  return (data.psd)
}
