# psdR
Power Spectral Density (psd) preprocessing method in R

### Usage

For now, the repo is private - to use the code
```
git clone git@github.com:abhivij/psdR.git
cd psdR
```

In R,
```
devtools::install()
library(psdR)
```

3 functions made available through this package :
* psd
* complexity
* compare_methods

Within R, to get man page for a function, say `complexity`
```
?complexity
```

The function definitions are present within the R directory.


The file `test.R`, contains sample function calls. It also contains, as comments, the execution time comparisons of `psd` in R and through Python.
These will be removed later.


### Usage in future

Once the repo is made public, the package can be installed by
```
devtools::install_github("abhivij/psdR")
```
