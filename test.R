# devtools::create("~/UNSW/VafaeeLab/psdR")

setwd("~/UNSW/VafaeeLab/psdR")
library(SingleCellExperiment)
library(reticulate)
source_python('psd.py')


# devtools::document()
#
# devtools::install()

compare_psd_implementation <- function(data, show_result = FALSE){
  start <- Sys.time()
  psd_r_data <- psdR::psd(data)
  end <- Sys.time()
  print(paste('Using R, psd execution time :', difftime(end, start, units = 'secs')))
  if (show_result) {
    print('psd output :')
    print(psd_r_data)
  }

  start <- Sys.time()
  psd_py_data <- psd(data)
  end <- Sys.time()
  print(paste('Using python via reticulate, psd execution time :', difftime(end, start, units = 'secs')))
  if (show_result) {
    print('psd output :')
    print(psd_py_data)
  }
  dimnames(psd_py_data) <- dimnames(psd_r_data)
  print(all.equal(psd_r_data, psd_py_data))
}

data1 <- matrix(rnorm(10), nrow=5, ncol=2)
compare_psd_implementation(data1)
# [1] "Using R, psd execution time : 0.000609159469604492"
# [1] "Using python via reticulate, psd execution time : 0.00504446029663086"
# [1] TRUE

data2 <- matrix(rnorm(9), nrow=3, ncol=3)
compare_psd_implementation(data2)
# [1] "Using R, psd execution time : 0.000649929046630859"
# [1] "Using python via reticulate, psd execution time : 0.00477814674377441"
# [1] TRUE

data <- matrix(c(-1.2225888, 1.34042501, 0.4259692,-0.23716855,
                  1.7332330, 0.21851021, -1.0023597, 1.00827522,
                  -0.6907052, 0.06109564),
                nrow=5, ncol=2, byrow = TRUE)
compare_psd_implementation(data)
# [1] "Using R, psd execution time : 0.000904321670532227"
# [1] "Using python via reticulate, psd execution time : 0.00530862808227539"
# [1] TRUE


data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
x <- counts(data)
x <- x[rowSums(x) > 0, ]
x <- t(as.matrix(x))
compare_psd_implementation(x)
# [1] "Using R, psd execution time : 213.227695226669"
# [1] "Using python via reticulate, psd execution time : 192.042806386948"
# [1] TRUE
