## Resubmission

This is a resubmission. In this version I have:

* Removed the `Date` field from the DESCRIPTION file.

* Added reference to the R journal paper about **roahd** in 

  + the `Description` field of the DESCRIPTION file, 
  + a CITATION file for use with the `citation()` function,
  + the `Citation` section of the README file,
  + the vignette.

* Auto-generated README from R markdown.

* Corrected a bug in the visualization of amplitude outliers for multivariate
functional data in `fbplot()`.

## Test environments

* I use the Github Action script
[`check-standard`](https://github.com/r-lib/actions/blob/master/examples/check-standard.yaml)
provided by the RStudio team to setup automatic environment checks. This
workflow runs `R CMD check` via the **rcmdcheck** package on the three major OSs
(linux, macOS and Windows) with the current release version of R, and R-devel.

* I have checked the package locally on my machine:
  - macOS Catalina 10.15.5, R-release 4.0.1.

* I have checked the package on R-hub on the following platforms:
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  - Ubuntu Linux 16.04 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran

* I have checked on win-builder for R-devel.

## R CMD check results

### [`check-standard`](https://github.com/r-lib/actions/blob/master/examples/check-standard.yaml) GitHub Action results

* There were no ERRORs.
* There were no WARNINGs.
* There were no NOTEs.

### Local (macOS) results

* There were no ERRORs.

* There was 1 WARNING because I have not installed **qpdf** on my machine:

```
* checking data for ASCII and uncompressed saves ... OK
   WARNING
  ‘qpdf’ is needed for checks on size reduction of PDFs
```

* There were 1 NOTE:

```
* checking installed package size ... NOTE
    installed size is  5.1Mb
    sub-directories of 1Mb or more:
      data   2.9Mb
      doc    1.7Mb
```

### R-hub results

* There were no ERRORs.

* There were no WARNINGs.

* There was 3 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Aymeric Stamm <aymeric.stamm@math.cnrs.fr>’

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  Ieva (38:28)
  al (38:36)
  et (38:33)
  roahd (39:34)
* checking installed package size ... NOTE
  installed size is  5.1Mb
  sub-directories of 1Mb or more:
    data   2.9Mb
    doc    1.7Mb
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
              user system elapsed
mfData       1.508  0.817   8.049
cor_spearman 1.876  0.028   6.593
```

### win-builder results

* There were no ERRORs.
* There were no WARNINGs.
* There were 1 NOTE:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Aymeric Stamm <aymeric.stamm@math.cnrs.fr>'

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  Ieva (38:28)
  al (38:36)
  et (38:33)
  roahd (39:34)
```
