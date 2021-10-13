N <- 1e2
P <- 1e2
grid <- seq(0, 1, length.out = P)
Cov <- exp_cov_function(grid, alpha = 0.3, beta = 0.4)

Data <- list()
Data[[1]] <- generate_gauss_fdata(
  N,
  centerline = sin(2 * pi * grid),
  Cov = Cov
)
Data[[2]] <- generate_gauss_fdata(
  N,
  centerline = sin(2 * pi * grid),
  Cov = Cov
)
names <- paste0("id_", 1:nrow(Data[[1]]))
DG <- depthgram(Data, marginal_outliers = TRUE, ids = names)
plot(DG)
