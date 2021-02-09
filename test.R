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
