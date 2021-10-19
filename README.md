# psdR
Power Spectral Density (psd) preprocessing method in R

### Usage

```
devtools::install_github("VafaeeLab/psdR")
```

In R,
```
library(psdR)
```

3 functions made available through this package :
* psd
* complexity
* compare_methods

Within R, to get help page for a function, say `complexity`
```
?complexity
```

The function definitions are present within the `R` directory in this repo.


The file `test.R`, contains sample function calls. It also contains, as comments, the execution time comparisons of `psd` in R and through Python.

`tsne_plot.png` and `complexity_df.txt` contain sample results of `compare_methods`. 

