# Digital video platforms and phylogenetic modelling reveal parrot tool use is not rare

Using YouTube and phylogenetic modelling to identify tool use in parrots.

## Getting started

### Software requirements

This code has been tested on the following operating systems:

- Windows 10 Education 22H2 64-bit
- Linux 3.10.0-693.2.2.el7.x86_64 x86_64
- macOS Ventura 13.4

This code requires R (v4.2.1) and the following R package versions:

- ape (v5.6-2)
- brms (v2.16.0)
- cowplot (v1.1.1)
- dagitty (v0.3-1)
- future (v1.27.0)
- future.callr (v0.8.0)
- geiger (v2.0.10)
- ggdag (v0.2.6)
- ggridges (v0.5.3)
- ggstance (v0.3.6)
- magick (v2.7.3)
- papaja (v0.1.1.9001)
- phangorn (v2.9.0)
- phytools (v1.0-3)
- readxl (v1.4.1)
- reshape2 (v1.4.4)
- rethinking (v2.21)
- rstan (v2.26.13)
- scales (v1.2.1)
- sdamr (v0.2.0)
- tarchetypes (v0.7.1)
- targets (v0.13.5)
- tidybayes (v3.0.2)
- tidyverse (v1.3.2)

### Installation guide

To run this code, you will need to [install R](https://www.r-project.org/) and the following R packages:

```
install.packages(
    c("ape", "brms", "cowplot", "dagitty", "future",
      "future.callr", "geiger", "ggdag", "ggridges",
      "ggstance", "magick", "papaja", "phangorn", 
      "phytools", "readxl", "reshape2", "rstan", 
      "scales", "sdamr", "tarchetypes", "targets",
      "tidybayes", "tidyverse")
)
```

You will also require the [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) package, which can be installed using the following code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

Finally, you will require the [rethinking](https://github.com/rmcelreath/rethinking) package, which can be installed by following the instructions [here](https://github.com/rmcelreath/rethinking#installation) or using the following code:

```
install.packages("cmdstanr")
cmdstanr::install_cmdstan()
install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
devtools::install_github("rmcelreath/rethinking")
```

These installation steps should typically take no longer than half an hour on a normal desktop computer.

## Run demo

A short demo executes the following steps from the full analysis pipeline: (1) loads the data, (2) estimates phylogenetic signal for tool use across 100 posterior samples from the parrot phylogeny, and (3) combines these estimates into a single posterior distribution.

To run this short demo:

1. Download this code repository to your local machine using `git clone https://github.com/ScottClaessens/phyloParrot` or by downloading the .zip file from GitHub
2. Set the working directory to this code repository `setwd("myPath/phyloParrot")` on your local machine
3. Load the `targets` package with `library(targets)`
3. Run `tar_make(phySignal1)` in the command line
4. To load the model object, run `tar_load(phySignal1)` in the command line

This demo should typically take 5-10 minutes on a normal desktop computer.

## Run full analysis pipeline

The full analysis pipeline can be viewed using `targets::tar_visnetwork()`. Running this pipeline will reproduce all figures, tables, and quantitative results from the manuscript, as well as a PDF of the manuscript itself.

To run the full analysis pipeline:

1. Download this code repository to your local machine using `git clone https://github.com/ScottClaessens/phyloParrot` or by downloading the .zip file from GitHub
2. Set the working directory to this code repository `setwd("myPath/phyloParrot")` on your local machine
3. Load the `targets` package with `library(targets)`
4. To run all analyses, run `tar_make()` in the command line
5. To load individual targets into your environment, run `tar_load(targetName)` in the command line

Note that some of the models, especially for the cross-validation, will take a long time to run. The full analysis pipeline will take several days to run in series. You can run the pipeline in parallel using `tar_make_clustermq()` or `tar_make_future()` (see [here](https://books.ropensci.org/targets/hpc.html)).

## Help

Any issues, please email scott.claessens@gmail.com.

## Authors

Scott Claessens, scott.claessens@gmail.com
