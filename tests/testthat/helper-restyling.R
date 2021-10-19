withr::local_seed(1234)
N <- 100
P <- 100
grid <- seq(0, 1, length.out = P)
centerline <- sin(2 * pi * grid)
Cov <- exp_cov_function(grid, alpha = 0.2, beta = 0.3)
Data_restyling <- generate_gauss_fdata(N, centerline, Cov)
fD_restyling <- fData(grid, Data_restyling)
