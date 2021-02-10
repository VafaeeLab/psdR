setwd("~/UNSW/VafaeeLab")
library(SingleCellExperiment)

# devtools::create("psdR")

data <- readRDS(file = 'psdR/TabulaMuris_Heart_10X.rds')
# data
# counts(data)
# rowData(data)
# rownames(data)
#
# colData(data)
#

setwd("~/UNSW/VafaeeLab/psdR")
devtools::document()

devtools::install()

data1 <- matrix(rnorm(10), nrow=5, ncol=2)
psdR::psd(data1)

data2 <- matrix(rnorm(9), nrow=3, ncol=3)
psdR::psd(data2)

data3 <- matrix(c(-1.2225888, 1.34042501, 0.4259692,-0.23716855,
                  1.7332330, 0.21851021, -1.0023597, 1.00827522,
                  -0.6907052, 0.06109564),
                nrow=5, ncol=2, byrow = TRUE)

data4 <- matrix(c(-1.2225888, 1.34042501, 0.4259692,-0.23716855,
                  1.7332330, 0.21851021, -1.0023597, 1.00827522,
                  -0.6907052, 0.06109564),
                nrow=1, byrow = TRUE)
