# devtools::create("~/UNSW/VafaeeLab/psdR")

setwd("~/UNSW/VafaeeLab/psdR")
Sys.setenv(keep.source.pkgs = TRUE)
devtools::document()
devtools::install()
devtools::load_all()


library(SingleCellExperiment)
Sys.setenv(RETICULATE_PYTHON = "../.py_venv/bin/python")
library(reticulate)
source_python('psd.py')




compare_psd_implementation <- function(data, show_result = FALSE){
  # start <- Sys.time()
  # psd_r_data <- psdR::psd_fft_apply(data)
  # end <- Sys.time()
  # print(paste('Using R fft apply, psd execution time :', difftime(end, start, units = 'secs')))

  start <- Sys.time()
  psd_r_data <- psdR::psd(data)
  end <- Sys.time()
  print(paste('Using R, psd execution time :', difftime(end, start, units = 'secs')))

  if (show_result) {
    print('psd output :')
    print(psd_r_data)
  }

  start <- Sys.time()
  psd_py_data <- t(psd(t(data)))
  end <- Sys.time()
  print(paste('Using python via reticulate, psd execution time :', difftime(end, start, units = 'secs')))
  if (show_result) {
    print('psd output :')
    print(psd_py_data)
  }
  dimnames(psd_py_data) <- dimnames(psd_r_data)
  print(all.equal(psd_r_data, psd_py_data))
}

data1 <- matrix(rnorm(10), nrow=2, ncol=5)
compare_psd_implementation(data1)
# [1] "Using R fft apply, psd execution time : 0.00855898857116699"
# [1] "Using R, psd execution time : 0.000279664993286133"
# [1] "Using python via reticulate, psd execution time : 0.00247502326965332"
# [1] TRUE

data2 <- matrix(rnorm(9), nrow=3, ncol=3)
compare_psd_implementation(data2)
# [1] "Using R fft apply, psd execution time : 0.00578093528747559"
# [1] "Using R, psd execution time : 0.000316381454467773"
# [1] "Using python via reticulate, psd execution time : 0.00202488899230957"
# [1] TRUE

data <- matrix(c(-1.2225888, 1.34042501, 0.4259692,-0.23716855,
                  1.7332330, 0.21851021, -1.0023597, 1.00827522,
                  -0.6907052, 0.06109564),
                nrow=2, ncol=5, byrow = TRUE)
compare_psd_implementation(data)
# [1] "Using R fft apply, psd execution time : 0.0092620849609375"
# [1] "Using R, psd execution time : 0.00032496452331543"
# [1] "Using python via reticulate, psd execution time : 0.00273799896240234"
# [1] TRUE


data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
x <- counts(data)
x <- x[rowSums(x) != 0, ]
x <- as.matrix(x)
compare_psd_implementation(x)
# [1] "Using R, psd execution time : 210.639530658722"
# [1] "Using R sweep, psd execution time : 212.108340740204"
# [1] "Using python via reticulate, psd execution time : 192.281036138535"
# [1] TRUE

# [1] "Using R fft apply, psd execution time : 205.484230279922"
# [1] "Using R, psd execution time : 219.209401130676"

# [1] "Using R fft apply, psd execution time : 211.930328607559"
# [1] "Using R, psd execution time : 212.845811367035"
# [1] "Using python via reticulate, psd execution time : 191.988132476807"
# [1] TRUE

# [1] "Using R fft apply, psd execution time : 205.381190061569"
# [1] "Using R, psd execution time : 204.640120029449"
# [1] "Using python via reticulate, psd execution time : 191.78946185112"
# [1] TRUE

##########################################
# complexity

data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
class_colname <- "cell_type1"
x <- as.matrix(SingleCellExperiment::counts(data))
classes <- SingleCellExperiment::colData(data)[class_colname]
classes <- classes[colnames(x),]
start <- Sys.time()
x <- x[rowSums(x) != 0, ]
f1 <- complexity(x, classes)
end <- Sys.time()
print(paste('complexity score = ', f1))
print(paste('complexity execution time :', difftime(end, start, units = 'secs')))

complexity(x)
complexity(x, c())

###############################################


data <- readRDS(file = 'TabulaMuris_Heart_10X.rds')
class_colname <- "cell_type1"
x <- as.matrix(SingleCellExperiment::counts(data))
classes <- SingleCellExperiment::colData(data)[class_colname]
classes <- classes[colnames(x),]

complexity_df_1 <- compare_methods(x, classes,
                                       c('CPM', 'NormCPM', 'TMM', 'Linnorm', 'Seurat'))
complexity_df_2 <- compare_methods(x, classes,
                                       c('Seurat'))
complexity_df_3 <- compare_methods(x, classes,
                                 c('CPM'))
complexity_df_4 <- compare_methods(data = x, classes = classes, plot_file_name = "cpm_tsne.png")
complexity_df_5 <- compare_methods(data = x, classes = NA, plot_file_name = "noclass_cpm_tsne.png")
