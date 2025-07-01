# semla <img src="man/figures/logo.png" align="right" width="120"/>

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-1.3.2-blue.svg)](https://github.com/ludvigla/semla/releases) [![](https://img.shields.io/github/last-commit/ludvigla/semla.svg)](https://github.com/ludvigla/semla/commits/main) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit/) [![DOI](https://zenodo.org/badge/599058747.svg)](https://zenodo.org/badge/latestdoi/599058747) [![Publication - 10.1093/bioinformatics/btad626](https://img.shields.io/badge/Publication-10.1093%2Fbioinformatics%2Fbtad626-2ea44f)](https://doi.org/10.1093/bioinformatics/btad626)

<!-- badges: end -->

<br>

`semla` is an R package that collects useful tools for Spatially Resolved Transcriptomics data analysis and visualization.

If you are visiting our Github page, please find more information at our package [website](https://ludvigla.github.io/semla/). Here, you can find documentation of functions together with examples on how to use them, as well as tutorials showing how to use `semla` for analysis and visualization of 10x Visium data.

For a more scientific description of the package, check out our associated publication "<i>Semla: a versatile toolkit for spatially resolved transcriptomics analysis and visualization</i>" in Bioinformatics ([DOI:10.1093/bioinformatics/btad626](https://doi.org/10.1093/bioinformatics/btad626)). Please cite this article if you use `semla` in your studies, and we'd be very grateful!

<br>

> [!NOTE]\
> Regarding [Seurat v5](https://satijalab.org/seurat/), `semla` functionalities should work as expected. However, if you run into any problems, please open an issue.

<br>

## Installation

The dev version of the package can be installed through GitHub using;

```         
install.packages("remotes")
remotes::install_github("ludvigla/semla")
```

<br>

### Setting up a conda environment

If you are using a [conda](https://docs.conda.io/en/latest/miniconda.html), you can follow the steps below:

```         
conda create -n r-semla r-essentials r-base
```

Activate the environment:

```         
conda activate r-semla
```

Then you can open RStudio from the environment. This should make sure that RStudio uses the R version and packages that are located in the conda environment. On Mac OS, you can open RStudio by typing the path to the RStudio application in the terminal:

```         
/Applications/RStudio.app/Contents/MacOS/RStudio
```

Note that `semla` requires R v4.1 or higher, so you should double check that you have the correct version installed.

Now you should be able to install `semla` in your conda environment:

```         
install.packages("remotes")
remotes::install_github("ludvigla/semla")
```

The `magick` R package might fail to install in your conda environment and in that case you can install it with conda instead:

```         
conda install -c conda-forge r-magick
```

<br>

### Renv

If you are familiar with [renv](https://rstudio.github.io/renv/articles/renv.html), you can install all necessary R packages (with exact versions) using `renv::restore()` with the `renv.lock` file provided in our GitHub repo.

<br>

### Docker

As an additional option, we provide a docker image on [Docker Hub](https://hub.docker.com/r/ludlar/semla) that you can download to run a container from. The image is based on the [rocker](https://hub.docker.com/r/rocker/rstudio) RStudio Server image and comes with RStudio Server pre-installed.

To access the server we need to publish a port for our container, we’ll use the --publish flag (-p for short) on the docker run command. The format of the --publish command is [host port]:[container port]. So, if we wanted to expose port 8000 inside the container to port 8080 outside the container, we would pass 8080:8000 to the --publish flag. We can also expose multiple ports which will be useful to use some of the interactive features in `semla` (see our `Interactive Viewer` section [here](https://ludvigla.github.io/semla/articles/feature_viewer.html)).

`<YOURPASSWORD>` can be any password you want to use for RStudio Server and the default username is `rstudio`. `<SOURCEPATH>` should be a path on your local machine which will be accessible from the docker container at `/home/rstudio`. Use `source="$(pwd)"` for current working directory.

```         
sudo docker run -d -p 1337:8787 -p 3030:3030 --name semla -e PASSWORD=<YOURPASSWORD> --memory=16g --mount type=bind,source="<SOURCEPATH>",target=/home/rstudio -e ROOT=TRUE ludlar/semla:latest
```

Visit the docker [docs](https://docs.docker.com/language/java/run-containers/) for more information on how to run containers.

<br>

## Citing *semla*

Larsson L, Franzén L, Ståhl PL, Lundeberg J. Semla: a versatile toolkit for spatially resolved transcriptomics analysis and visualization. Bioinformatics. 2023 Oct 3;39(10):btad626. doi: 10.1093/bioinformatics/btad626. PMID: 37846051; PMCID: PMC10597621.
https://academic.oup.com/bioinformatics/article/39/10/btad626/7319366 

<br>

## What is semla?

A semla is a delicious pastry, traditionally consumed around a specific day, "Semmeldagen" or "Fettisdagen" ("Fat Tuesday"), in Sweden. It is made up of a light cardamom flavored wheat bun, fluffy whipped cream, sweet crunchy almond paste, and topped with dusted icing sugar. If that's not reason enough for us to name this package after this fantastic creation, let's also say *semla* is an abbreviation for "**S**patially r**E**solved transcripto**M**ics too**L**s for **A**nalysis".

<br>
