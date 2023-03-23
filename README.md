# semla <img src="vignettes/logo.png" align="right" height="120"/>

<!-- badges: start -->

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![](https://img.shields.io/github/last-commit/ludvigla/semla.svg)](https://github.com/ludvigla/semla/commits/main)

<!-- badges: end -->

`semla` is an R package that collects useful tools for Spatially Resolved Transcriptomics data analysis and visualization.

If you are visiting our Github page, please find more information at our package [website](https://ludvigla.github.io/semla/)

Here you can find documentation of functions together with examples on how to use them, as well as tutorials showing how to use `semla` for analysis of 10x Visium data.

## Installation

The dev version of the package can be installed through GitHub using;

````
remotes::install_github("ludvigla/semla")
````


## Setting up a conda environment

If you want to run `semla` within a conda environment, you first need to ensure to have anaconda installed. I suggest using [miniconda](https://docs.conda.io/en/latest/miniconda.html)

    conda create -n r-semla r-essentials r-base

There's one R package that might cause a few issues if not configured properly, namely `magick`. I suggest installing it with conda within the environment:

    conda install -c conda-forge r-magick

Now activate the environment

    conda activate r-semla

Then you can open RStudio from the environment. This should make sure that RStudio uses the R version and packages that are located in the conda environment. On Mac OS, you can open RStudio by running something like:

    /Applications/RStudio.app/Contents/MacOS/RStudio

When RStudio is opened, check that the R version is higher than v4.1.

## What is semla?

A semla is a delicious pastry, traditionally consumed around a specific day, "Semmeldagen" or "Fettisdagen" ("Fat Tuesday"), in Sweden. It is made up of a light cardamom flavored wheat bun, fluffy whipped cream, sweet crunchy almond paste, and topped with dusted icing sugar. If that's not reason enough for us to name this package after this fantastic creation, let's also say *semla* is an abbreviation for "**S**patially r**E**solved transcripto**M**ics too**L**s for **A**nalysis".

<br>
