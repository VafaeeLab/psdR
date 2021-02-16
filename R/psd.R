#' psd function
#'
#' Compute the Power Spectral Density (psd) transformation of the given scRNA count data
#'
#' @param data scRNA count data - cells Vs genes (dim : c x g)
#' @return psd transformation of input data
#' @export
psd <- function(data){
  data.psd <- cor(data)                               #dim : g x g
  data.psd <- data %*% data.psd                       #dim : c x g
  data.psd <- abs(t(                                  #compute discrete fourier transform
                    mvfft(t(abs(data.psd)))                         #of each cell sample
                   )
                 ) / dim(data)[2]
  data.psd <- data.psd / rowSums(data.psd)            #divide rowSums from each element
  data.psd <- data.psd * log2(data.psd)
  data.psd <- apply(data.psd, 2,
                    function(x){return (x - mean(x))}
                    )                                 #subtract column mean from each element
  data.psd <- apply(data.psd, 2,
                    function(x){
                      return ( (x-min(x)) / (max(x)-min(x)) )
                    })
  return (data.psd)
}
