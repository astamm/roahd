# Changelog



Here's a list of what is changed in this update of __roahd__:

### Major fixes 

1) Removed check for uniformity in the grid of fData() and mfData() constructor

2) Added the possibility to subset fData in time with logical vectors

3) Fixes in methods BD, BD_relative, HI and EI: the previous computational technique was based on arguments from the popular reference "Exact fast computation of band depth for large functional datasets: How quickly can one million curves be ranked?" by Sun, Genton and Nychka, which in the case of BD, and HI/EI are wrong. Now the implementation exploited sticks to the definition, at the cost of a higher computational burden (and thus, time to complete the computation).
