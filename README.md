# pseudogp

pseudogp is an R package for Bayesian inference of Gaussian Process Latent Variable models learning pseudotimes from single-cell RNA-seq. It forms a wrapper round a [stan](http://mc-stan.org/) model (that requires the [rstan package](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)), as well as a set of functions to plot the posterior mean curves.

pseudogp acts on a reduced dimension representation of the data (in a similar manner to [monocle](http://cole-trapnell-lab.github.io/monocle-release/) and [waterfall](http://www.cell.com/cell-stem-cell/fulltext/S1934-5909(15)00312-4)) in which it fits a probabilsitic curve, allowing posterior pseudotime uncertainty to be quantified.



## Installation

```R
# install.packages("devtools")
devtools::install_github("kieranrcampbell/pseudogp")
```
