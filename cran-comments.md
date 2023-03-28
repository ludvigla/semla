## CRAN submission

Please omit the redundant "Tools for" in your title and "A toolkit for" from the start of your description text.

__Comment: Done.__

Unexecutable code in man/manual-transform-images.Rd
-> Please check that

__Comment: Fixed example.__

Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
Missing Rd-tags in up to 11 .Rd files, e.g.:
export_graph.Rd: \value
ftrviewer-shiny.Rd: \value
ftrviewer.Rd: \value
osddu-shiny.Rd: \value
osddu.Rd: \value
paper-shiny.Rd: \value

__Comment: Added \value to exported functions in .Rd files.__

Some code lines in examples are commented out.
Please never do that. Ideally find toy examples that can be regularly executed and checked. Lengthy examples (> 5 sec), can be wrapped in \donttest{}.
Examples in comments in:
UpdateSTUtilityV1Object.Rd

__Comment: Removed comment from example in UpdateSTUtilityV1Object.Rd.__

\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}.

__Comment: Removed dontrun from some tests and replaces with donttest for the remaining ones.__

Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies.
Please omit any default path in writing functions. In your examples/vignettes/tests you can write to tempdir().

__Comment: Removed default paths from functions and changed to tempdir() in examples.__

Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)  # code line i
on.exit(par(oldpar))  # code line i + 1
...
par(mfrow=c(2,2))  # somewhere after
...
e.g.: R/plot_images.R
If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.

Please always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos. -> man/load-images.Rd ; man/mask-images.Rde.g.:
oldpar <- par(mfrow = c(1,2))

__Comment: Changed code in R/plot_images.R to return to previous graphical parameter settings on exit. Also changed examples in load_images.R and mask.R.__

Please do not install packages in your functions, examples or vignette. This can make the functions,examples and cran-check very slow.

__Comment: Removed calls to install.packages from examples and functions.__

Please ensure that you do not use more than 2 cores in your examples, vignettes, etc.

__Comment: TBA.__


## CRAN pre-tests (2023-03-27)

Found the following (possibly) invalid URLs:
  URL: https://opensource.org/licenses/mit/ (moved to https://opensource.org/license/mit/)
    From: README.md
    Status: 301
    Message: Moved Permanently
    
Comment: Fixed typo in URL.

## CRAN pre-tests

* checking CRAN incoming feasibility ... [9s] NOTE
Possibly misspelled words in DESCRIPTION:
  Transcriptomics (3:37, 24:47)
  semla (24:94)
  
Comment: Transcriptomics is not misspelled. Added single quotes to package name semla.

Found the following (possibly) invalid URLs:
  URL: https://opensource.org/licenses/mit (moved to https://opensource.org/license/mit/)
    From: README.md
    Status: 301
    Message: Moved Permanently
    
Comment: Added trailing slash to URL.

* checking package dependencies ... NOTE
Package suggested but not available for checking: 'TabulaMurisSenisData'

Comment: package not hosted on CRAN

Size of tarball: 6165189 bytes

Not more than 5MB for a CRAN package.

Comment: Reduced size of example data. Tarball size is now 3747730 bytes.


## Test environments

* local macOS Big Sur 10.16 install (devel), R version 4.2.1
* local macOS Catalina 10.15.7 install (devel), R version 4.1.2
* local macOS Monterey 12.3.1 install (devel), R version 4.2.1
* local Windows 10 Education, 64-bit
* rhub (Fedora Linux, R-devel, clang, gfortran)
* rhub (Ubuntu Linux 20.04.1 LTS, R-release, GCC)
* rhub (Windows Server 2022, R-devel, 64 bit)

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 1 note ✖

1 NOTE:
❯ checking installed package size ... NOTE
    installed size is  8.3Mb
    sub-directories of 1Mb or more:
      extdata       5.6Mb
      htmlwidgets   1.5Mb

Comment: 
extdata/ contains two example datasets which have been reduced to a fraction of the 
original data size. These example datasets are required for the function examples and 
the tutorials provided on our pkgdown [website](https://github.com/ludvigla/semla) and
we hope that the data is sufficiently small to stay in the package as it is.
htmlwidgets/ contains javascript code required for interactive web-applications provided 
in the package and cannot be removed.

## Reverse dependencies

None