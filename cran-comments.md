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

* There were no NOTEs.

### R-hub results

* There were no ERRORs.

* There were no WARNINGs.

* There was 1 NOTE because the package was previously on CRAN but has been
removed automatically due to incompatibility with new R-release 4.0 and because
of a change in maintainer:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Aymeric Stamm <aymeric.stamm@math.cnrs.fr>’

New submission

Package was archived on CRAN

The Date field is over a month old.
```

### win-builder results

* There were no ERRORs.
* There were no WARNINGs.
* There were no NOTEs.
