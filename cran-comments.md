## Test environments
* local macOS R installation, R 4.1.2
* macOS latest release (via [R-CMD-check](https://github.com/r-lib/actions/blob/master/examples/check-standard.yaml) github action)
* windows latest release (via [R-CMD-check](https://github.com/r-lib/actions/blob/master/examples/check-standard.yaml) github action)
* ubuntu 20.04 latest both release and devel (via [R-CMD-check](https://github.com/r-lib/actions/blob/master/examples/check-standard.yaml) github action)
* [win-builder](https://win-builder.r-project.org/) (release and devel)
* [R-hub](https://builder.r-hub.io)
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

    * checking installed package size ... NOTE
      installed size is  5.2Mb
      sub-directories of 1Mb or more:
        data   2.9Mb
        doc    1.7Mb

The size varies according to the system on which the package is installed.

    * Maintainer: 'Aymeric Stamm <aymeric.stamm@math.cnrs.fr>'
      Possibly misspelled words in DESCRIPTION:
        Aleman (42:55)
        depthgram (43:57)

This is due to a change of maintainer from Nicholas Tarabelloni to Aymeric Stamm.
Words listed as possibly misspelled are people name and name of a proposed mathematical method, hence they are spelled correctly.
