# Changelog

Here's a list of what is changed in this update of __roahd__:

### Upgrades

#### Major upgrades

1) Extended Spearman's correlation coefficient computation for multivariate datasets with more than two 
components.

2) Added bootstrap-based computation of Spearman's correlation coefficient bias and standard deviation.

3) Added methods to provide bootstrap-based confidence intervals on Spearman's coefficients for two 
univariate functional datasets or a multivairate functional dataset.

4) Added a bootstrap-based test on Spearman's correlation coefficient for two multivariate functional datasets.

5) Added an outliergram version (without graphical display of original data) of multivariate functional datasets.

6) Added example multivariate functional datasets of ECG signals.

#### Minor updates

1) Added two convenience functions to append compatible functional datasets (univariate or multivariate).

2) Added a [-operator overload for multivariate functional dataset representation __mfData__.


### Fixes

#### Major fixes

1) Fixed bug in cor_spearman function. Now the standard spearman correlation is not computed on ranks of MHI/MEI, but on
MHI/MEI itself. The difference is very small, but allows for full reproducibility of the results in the original paper.

#### Minor fixes 

1) Fixed typos in doc

2) Standardised formulas for the application of F inflations in outliergram and boxplot

