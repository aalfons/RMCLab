# rdmc: Robust Discrete Matrix Completion

Robust matrix completion for discrete rating-scale data with a low-rank constraint on a latent continuous matrix, implemented via an ADMM algorithm. For the loss part of the objective function, several robust loss functions are implemented. In addition, the package provides wrapper functions for `softImpute` (Mazumder, Hastie, and Tibshirani, 2010, <https://www.jmlr.org/papers/v11/mazumder10a.html>; Hastie, Mazumder, Lee, Zadeh, 2015, <https://www.jmlr.org/papers/v16/hastie15a.html>) for easy tuning of the regularization parameter, as well as benchmark methods such as median imputation and mode imputation.


## Installation

To install the latest version from GitHub, you can pull this repository and 
install it from the `R` command line via

```
install.packages("devtools")
devtools::install_github("aalfons/rdmc")
```

If you already have package `devtools` installed, you can skip the first
line. Moreover, package `rdmc` contains `C++` code that needs to be
compiled, so you may need to download and install the [necessary tools
for MacOS](https://cran.r-project.org/bin/macosx/tools/) or the
[necessary tools for
Windows](https://cran.r-project.org/bin/windows/Rtools/).


## Report issues and request features

If you experience any bugs or issues or if you have any suggestions for additional features, please submit an issue via the [*Issues*](https://github.com/aalfons/rdmc/issues) tab of this repository.  Please have a look at existing issues first to see if your problem or feature request has already been discussed.


## Contribute to the package

If you want to contribute to the package, you can fork this repository and create a pull request after implementing the desired functionality.


## Ask for help

If you need help using the package, or if you are interested in collaborations related to this project, please get in touch with the [package maintainer](https://personal.eur.nl/alfons/).
