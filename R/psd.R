#' psd function
#'
#' Compute the Power Spectral Density (psd) transformation of the given scRNA count data
#'
#' @param data scRNA count data - genes Vs cells (dim : g x c)
#' @return psd transformation of input data (dim : g x c)
#' @export
psd <- function(data){
  data <- t(data)                           #dim : c x g
  data.psd <- cor(data)                     #dim : g x g
  data.psd <- data %*% data.psd             #dim : c x g
  data.psd <- abs(mvfft(t(abs(data.psd)))   #compute discrete fourier transform of each cell sample
                 ) / dim(data)[2]           #dim : g x c
  data.psd <- apply(data.psd, 2,            #divide column sum from each element
                    function(x){x / sum(x)}
                    )
  data.psd <- data.psd * log2(data.psd)
  data.psd <- apply(data.psd, 1,            #subtract row mean from each element
                    function(x){return (x - mean(x))}
                    )
  data.psd <- apply(data.psd, 2,
                    function(x){
                      return ( (x-min(x)) / (max(x)-min(x)) )
                    })
  return (t(data.psd))
}


# data <- matrix(c(-1.2225888, 1.34042501, 0.4259692,-0.23716855,
#                  1.7332330, 0.21851021, -1.0023597, 1.00827522,
#                  -0.6907052, 0.06109564),
#                nrow=2, ncol=5, byrow = TRUE)
# psd(data)
