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