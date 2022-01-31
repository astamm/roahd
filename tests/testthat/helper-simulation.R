N <- 30
P <- 100
L <- 3

t0 <- 0
tP <- 1
time_grid <- seq(t0, tP, length.out = P)

C1 <- exp_cov_function(time_grid, alpha = 0.1, beta = 0.2)
C2 <- exp_cov_function(time_grid, alpha = 0.2, beta = 0.5)
C3 <- exp_cov_function(time_grid, alpha = 0.3, beta = 1)

CholC1 <- chol(C1)
CholC2 <- chol(C2)
CholC3 <- chol(C3)

centerline <- sin(2 * pi * time_grid)
centerlines <- matrix(c(
  centerline,
  sqrt(time_grid),
  10 * (time_grid - 0.5) * time_grid
  ), nrow = 3, byrow = TRUE)
