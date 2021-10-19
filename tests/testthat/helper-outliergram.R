withr::local_seed(1234)

N <- 200
N_outliers <- 4
P <- 200
grid <- seq(0, 1, length.out = P)
Cov = exp_cov_function(grid, alpha = 0.2, beta = 0.8)

Data <- generate_gauss_fdata(
  N,
  centerline = sin(4 * pi * grid),
  Cov = Cov
)

Data_out <- array(0, dim = c(N_outliers, P))
Data_out[1, ] <- generate_gauss_fdata(
  1,
  centerline = sin(4 * pi * grid + pi / 2),
  Cov = Cov
)
Data_out[2, ] <- generate_gauss_fdata(
  1,
  centerline = sin(4 * pi * grid - pi / 2),
  Cov = Cov
)
Data_out[3, ] <- generate_gauss_fdata(
  1,
  centerline = sin(4 * pi * grid + pi/ 3),
  Cov = Cov
)
Data_out[4, ] <- generate_gauss_fdata(
  1,
  centerline = sin(4 * pi * grid - pi / 3),
  Cov = Cov
)

Data <- rbind(Data, Data_out)
fDout <- fData(grid, Data)
