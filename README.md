# Crowdsourcing and phylogenetic modelling reveal parrot tool use is not rare

Using crowdsourcing and phylogenetic modelling to identify tool use in parrots.

## Getting Started

### Downloading data

All data for this project can be found in the `data/` folder. To download the data, either `git clone` this repository or download the repository as a .zip file.

### Installing

To run this code, you will need to [install R](https://www.r-project.org/) and the following R packages:

```
install.packages(
    c("ape", "brms", "cowplot", "dagitty", "geiger",
      "ggdag", "ggridges", "ggstance", "magick",
      "papaja", "phangorn", "phytools", "readxl",
      "reshape2", "rethinking", "rstan", "scales",
      "sdamr", "tidybayes", "tidyverse")
)
```

You will also require the [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) package, which can be installed using the following code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

### Executing code

1. Set the working directory to this code repository `setwd("myPath/phyloParrot")` on your local machine
2. Load the `targets` package with `library(targets)`
3. To run all analyses, run `tar_make()`
4. To load individual targets into your environment, run `tar_load(targetName)`

Note that some of the models, especially for the cross-validation, will take a long time to run. You can run the pipeline in parallel using `tar_make_clustermq()` or `tar_make_future()` (see [here](https://books.ropensci.org/targets/hpc.html)).

## Help

Any issues, please email scott.claessens@gmail.com.

## Authors

Scott Claessens, scott.claessens@gmail.com
