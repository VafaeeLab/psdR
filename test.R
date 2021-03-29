# devtools::create("~/UNSW/VafaeeLab/psdR")

setwd("~/UNSW/VafaeeLab/psdR")
devtools::document()
devtools::install()


library(SingleCellExperiment)
Sys.setenv(RETICULATE_PYTHON = "../.py_venv/bin/python")
library(reticulate)
source_python('psd.py')




compare_psd_implementation <- function(data, show_result = FALSE){
  start <- Sys.time()
  psd_r_data <- psdR::psd(data)
  end <- Sys.time()
  print(paste('Using R, psd execution time :', difftime(end, start, units = 'secs')))

  start <- Sys.time()
  psd_r_data <- psdR::psd_without_t(data)
  end <- Sys.time()
  print(paste('Using R sweep, psd execution time :', difftime(end, start, units = 'secs')))

  if (show_result) {
    print('psd output :')
    print(psd_r_data)
  }

  start <- Sys.time()
  psd_py_data <- psd(t(data))
  end <- Sys.time()
  print(paste('Using python via reticulate, psd execution time :', difftime(end, start, units = 'secs')))
  psd_py_data <- t(psd_py_data)
  if (show_result) {
    print('psd output :')
    print(psd_py_data)
  }
  dimnames(psd_py_data) <- dimnames(psd_r_data)
  print(all.equal(psd_r_data, psd_py_data))
}

data1 <- matrix(rnorm(10), nrow=2, ncol=5)
compare_psd_implementation(data1)
# [1] "Using R, psd execution time : 0.00162935256958008"
# [1] "Using python via reticulate, psd execution time : 0.00282073020935059"
# [1] TRUE

data2 <- matrix(rnorm(9), nrow=3, ncol=3)
compare_psd_implementation(data2)
# [1] "Using R, psd execution time : 0.00117683410644531"
# [1] "Using python via reticulate, psd execution time : 0.0077660083770752"
# [1] TRUE

data <- matrix(c(-1.2225888, 1.34042501, 0.4259692,-0.23716855,
                  1.7332330, 0.21851021, -1.0023597, 1.00827522,
                  -0.6907052, 0.06109564),
                nrow=2, ncol=5, byrow = TRUE)
compare_psd_implementation(data)
# [1] "Using R, psd execution time : 0.00132465362548828"
# [1] "Using python via reticulate, psd execution time : 0.00522065162658691"
# [1] TRUE


data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
x <- counts(data)
x <- x[rowSums(x) > 0, ]
x <- as.matrix(x)
compare_psd_implementation(x)
# [1] "Using R, psd execution time : 210.639530658722"
# [1] "Using R sweep, psd execution time : 212.108340740204"
# [1] "Using python via reticulate, psd execution time : 192.281036138535"
# [1] TRUE



##########################################
# complexity

data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
class_colname <- "cell_type1"
start <- Sys.time()
f1 <- psdR::complexity(data, class_colname)
end <- Sys.time()
print(paste('complexity score = ', f1))
print(paste('complexity execution time :', difftime(end, start, units = 'secs')))
