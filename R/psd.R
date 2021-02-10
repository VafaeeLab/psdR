#' psd function
#'
#' Compute the Power Spectral Density (psd) transformation of the given scRNA count data
#'
#' @param data scRNA count data - cells Vs genes
#' @return psd transformation of input data
#' @export
psd <- function(data){
  data.cor <- cor(data)               #dim : g x g
  data.psd <- data %*% data.cor       #dim : c x g
                #compute discrete fourier transform of each cell sample
  data.psd <- abs(t(
                  mvfft(t(abs(data.psd)))
                  )) / dim(data)[2]

  data.psd <- data.psd / rowSums(data.psd)            #subtract rowSums from each column
  data.psd <- data.psd * log2(data.psd)
  data.psd <- sweep(data.psd, 2, colMeans(data.psd))  #subtract colMeans from each row

  # below line doesn't give same results as in psd.py ! Why ?
  # data.psd <- apply(data.psd, 1, min_max_normalize)

  col_min <- apply(data.psd, 2, min)
  col_max <- apply(data.psd, 2, max)

  data.psd <- sweep(sweep(data.psd, 2, col_min), 2, col_max - col_min, "/")

  return (data.psd)
}

# min_max_normalize <- function(x){
#   return (x - min(x)) / (max(x) - min(x))
# }
