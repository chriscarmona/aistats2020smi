
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Semi-Modular Inference

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![devel-version](https://img.shields.io/badge/devel%20version-0.2.0-blue.svg)](https://github.com/christianu7/aistats2020smi)
<!-- badges: end -->

This repo contains code for our AISTATS article on **Semi-Modular
Inference**.

Semi-Modular Inference (SMI) is a modification of Bayesian inference in
multi-modular settings, which enables tunable and directed flow of
information between modules.

For an introduction to SMI, we invite you to watch our [slideslive
presentation](https://slideslive.com/38930337) (best on 1.5x),
[<img src="inst/figures/cover.png" width="70%" style="display: block; margin: auto;" />](https://slideslive.com/38930337)

## Citation

If you find Semi-Modular Inference relevant for your scientific
publication, we encourage you to add the following reference:

``` bibtex
@InProceedings{Carmona2020smi,
  title = {Semi-Modular Inference: enhanced learning in multi-modular models by tempering the influence of components},
  author = {Carmona, Chris U. and Nicholls, Geoff K.},
  booktitle = {Proceedings of the 23rd International Conference on Artificial Intelligence and Statistics, AISTATS 2020},
  year = {2020},
  editor = {Silvia Chiappa and Roberto Calandra},
  volume = {108},
  pages = {4226--4235},
  series = {Proceedings of Machine Learning Research},
  month = {26--28 Aug},
  publisher = {PMLR},
  pdf = {http://proceedings.mlr.press/v108/carmona20a/carmona20a.pdf},
  url = {http://proceedings.mlr.press/v108/carmona20a.html},
  arxivId = {2003.06804},
}
```

## Installation

You can install the devel version of aistats2020smi from our github
repository

``` r
install_github("christianu7/aistats2020smi")
```

## Reproducibility

The main article and supplementary material can be reproduced entirely
using a `.Rnw` file included in this repo. Executing the following
command will generate a pdf file in your `Documents` folder:

``` r
aistats2020smi::generate_article( out_dir="~/Documents" )
```

If you prefer to keep and analyse intermediate outputs, consider
executing the following commands:

``` r
dir.create("~/smi_article")
aistats2020smi::download_mcmc_results( dest_path = "~/smi_article" )
aistats2020smi::generate_article( out_dir="~/smi_article", mcmc_dir="~/smi_article" )
```
