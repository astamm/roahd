#' Depthgram for univariate and multivariate functional data sets
#'
#' This function computes the three 'DepthGram' representations from a p-variate
#' functional data set.
#'
#' @param Data A \code{\link[base]{list}} of length `L` (number of components)
#'   in which each element is an `N x P` matrix with `N` individuals and `P`
#'   time points. Alternatively, it can also be an object of class
#'   \code{\link{fData}} or of class \code{\link{mfData}}.
#' @param marginal_outliers A boolean specifying whether the function should
#'   return shape and amplitude outliers over each dimension. Defaults to
#'   `FALSE`.
#' @param boxplot_factor A numeric value specifying the inflation factor for
#'   marginal functional boxplots. This is ignored if `marginal_outliers ==
#'   FALSE`. Defaults to `1.5`.
#' @param outliergram_factor A numeric value specifying the inflation factor for
#'   marginal outliergrams. This is ignored if `marginal_outliers == FALSE`.
#'   Defaults to `1.5`.
#' @param ids A character vector specifying labels for individual observations.
#'   Defaults to `NULL`, in which case observations will remain unlabelled.
#'
#' @return An object of class `depthgram` which is a list with the following
#'   items:
#'
#' - `mbd.mei.d`: vector MBD of the MEI dimension-wise.
#' - `mei.mbd.d`: vector MEI of the MBD dimension-wise.
#' - `mbd.mei.t`: vector MBD of the MEI time-wise.
#' - `mei.mbd.t`: vector MEI of the MEI time-wise.
#' - `mbd.mei.t2`: vector MBD of the MEI time/correlation-wise.
#' - `mei.mbd.t2`: vector MEI of the MBD time/correlation-wise.
#' - `shp.out.det`: detected shape outliers by dimension.
#' - `mag.out.det`: detected magnitude outliers by dimension.
#' - `mbd.d`: matrix `n x p` of MBD dimension-wise.
#' - `mei.d`: matrix `n x p` of MEI dimension-wise.
#' - `mbd.t`: matrix `n x p` of MBD time-wise.
#' - `mei.t`: matrix `n x p` of MEI time-wise.
#' - `mbd.t2`: matrix `n x p` of MBD time/correlation-wise
#' - `mei.t2`: matrix `n x p` of MBD time/correlation-wise.
#'
#' @references
#' Aleman-Gomez, Y., Arribas-Gil, A., Desco, M. Elias-Fernandez, A., and Romo,
#' J. (2021). "Depthgram: Visualizing Outliers in High Dimensional Functional
#' Data with application to Task fMRI data exploration".
#'
#' @export
#'
#' @examples
#' N <- 2e2
#' P <- 1e3
#' grid <- seq(0, 1, length.out = P)
#' Cov <- exp_cov_function(grid, alpha = 0.3, beta = 0.4)
#'
#' Data <- list()
#' Data[[1]] <- generate_gauss_fdata(
#'   N,
#'   centerline = sin(2 * pi * grid),
#'   Cov = Cov
#' )
#' Data[[2]] <- generate_gauss_fdata(
#'   N,
#'   centerline = sin(2 * pi * grid),
#'   Cov = Cov
#' )
#' names <- paste0("id_", 1:nrow(Data[[1]]))
#'
#' DG1 <- depthgram(Data, marginal_outliers = TRUE, ids = names)
#'
#' fD <- fData(grid, Data[[1]])
#' DG2 <- depthgram(fD, marginal_outliers = TRUE, ids = names)
#'
#' mfD <- mfData(grid, Data)
#' DG3 <- depthgram(mfD, marginal_outliers = TRUE, ids = names)
depthgram <- function(Data,
                      marginal_outliers = FALSE,
                      boxplot_factor = 1.5,
                      outliergram_factor = 1.5,
                      ids = NULL) {
  UseMethod("depthgram", Data)
}

#' @rdname depthgram
#' @export
depthgram.default <- function(Data,
                              marginal_outliers = FALSE,
                              boxplot_factor = 1.5,
                              outliergram_factor = 1.5,
                              ids = NULL) {
  p <- length(Data)
  n <- nrow(Data[[1]])
  N <- ncol(Data[[1]])

  if (marginal_outliers) {
    a2 <- a0 <- -2 / (n * (n - 1))
    a1 <- 2 * (n + 1) / (n - 1)
  }

  #########################
  # Dimension-wise
  #########################

  # n*p matrix with mbd's on each dimension
  mbd.d <- array(0, dim = c(n, p))
  # n*p matrix with mbd's on each dimension
  mei.d <- array(0, dim = c(n, p))
  # list for marginal magnitude outliers detection
  mag.out.det <- list(length = p)
  # list for marginal shape outliers detection
  shp.out.det <- list(length = p)
  n2 <- ceiling(n * 0.5)

  # Vector containing the sign of correlation between mei in each dimension and
  # the next one
  corr.mei <- vector("numeric", length = p)
  corr.mei[1] <- 1

  #Array with all observation ranks on each time-point/dimension
  rmat.mat <- array(0, dim = c(n, N, p))
  # Storing components in which -1 transformation will be applied (since negative
  # mei correlation exists)
  wp <- c()

  for (i in 1:p) { ### Over dimensions
    x <- Data[[i]]
    # MBD and MEI computation on dimension i
    rmat <- apply(t(x), 1, rank)
    rmat.mat[, , i] <- rmat
    down <- rmat - 1
    up <- n - rmat
    mbd.d[, i] <- (rowSums(up * down) / N + n - 1) / (n * (n - 1) / 2)
    mei.d[, i] <- rowSums(up + 1) / (n * N)
    # MEI correlation between i and i-1 dimensions
    if (i > 1)
      corr.mei[i] <- corr.mei[i - 1] * sign(stats::cor(mei.d[, i], mei.d[, i - 1]))

    if (corr.mei[i] == -1)
      wp <- c(wp, i)

    # Marginal outlier detection on dimension i: functional boxplot and
    # outliergram with factors boxplot_factor and outliergram_factor
    # respectively
    if (marginal_outliers) {
      index <- order(mbd.d[, i], decreasing = TRUE)
      center <- x[index[1:n2], ]
      inf <- apply(center, 2, min)
      sup <- apply(center, 2, max)
      dist <- boxplot_factor * (sup - inf)
      upper <- sup + dist
      lower <- inf - dist
      mag.out.det[[i]] <- which(colSums((t(x) <= lower) + (t(x) >= upper)) > 0)
      dist <- (a0 + a1 * mei.d[, i] + a2 * n^2 * mei.d[, i]^2) - mbd.d[, i]
      q <- stats::quantile(dist, probs = c(0.25, 0.75))
      lim <- outliergram_factor * (q[2] - q[1]) + q[2]
      shp.out.det[[i]] <- which(dist > lim)
      rm(index, center, inf, sup, dist, upper, lower, q, lim)
    }

    rm(x, rmat, down, up)
  }

  rm(Data)

  ##### MEI of MBD and MBD of MEI across dimensions

  mei.mbd.d <- MEI(mbd.d) # managing ties
  mbd.mei.d <- MBD(mei.d)

  ########################
  # Time-wise
  ########################

  # n*N matrix with mbd's on each time point
  mbd.t  <- array(0, dim = c(n, N))
  # n*N matrix with mei's on each time point
  mei.t  <- array(0, dim = c(n, N))
  # n*N matrix with mbd's on each time point for the "corrected" data set
  mbd.t2 <- array(0, dim = c(n, N))
  # n*N matrix with mbd's on each time point for the "corrected" data set
  mei.t2 <- array(0, dim = c(n, N))

  for (i in 1:N) { ### Over time points
    # Getting observation ranks at time-point i
    rmat <- rmat.mat[, i, ]
    if (is.null(dim(rmat)))
      rmat <- matrix(rmat, ncol = 1)
    # MBD and MEI computation on time-point i
    down <- rmat - 1
    up <- n - rmat
    mbd.t[, i] <- (rowSums(up * down) / p + n - 1) / (n * (n - 1) / 2)
    mei.t[, i] <- rowSums(up + 1) / (n * p)

    # MBD and MEI computation on dimension i for the "corrected" data set
    if (length(wp) > 0) {
      down[, wp] <- n - rmat[, wp]
      up[, wp] <- rmat[, wp] - 1
      mbd.t2[, i] <- (rowSums(up * down) / p + n - 1) / (n * (n - 1) / 2)
      mei.t2[, i] <- rowSums(up + 1) / (n * p)
    } else {
      mbd.t2[, i] <- mbd.t[, i]
      mei.t2[, i] <- mei.t[, i]
    }

    rm(down, rmat, up)
  }

  ##### MEI of MBD and MBD of MEI across time-points (original and corrected
  ##### data sets)
  mei.mbd.t  <- MEI(mbd.t)
  mei.mbd.t2 <- MEI(mbd.t2)
  mbd.mei.t  <- MBD(mei.t)
  mbd.mei.t2 <- MBD(mei.t2)

  ##### RETURN
  res <- list(
    mbd.mei.d  = mbd.mei.d,  mei.mbd.d  = mei.mbd.d,
    mbd.mei.t  = mbd.mei.t,  mei.mbd.t  = mei.mbd.t,
    mbd.mei.t2 = mbd.mei.t2, mei.mbd.t2 = mei.mbd.t2,
    shp.out.det = shp.out.det, mag.out.det = mag.out.det,
    mbd.d  = mbd.d,  mei.d  = mei.d,
    mbd.t  = mbd.t,  mei.t  = mei.t,
    mbd.t2 = mbd.t2, mei.t2 = mei.t2
  )

  if (!is.null(ids)) {
    res <- lapply(res, function(.x) {
      if (is.list(.x)) {
        .x <- lapply(.x, function(.y) {
          names(.y) <- ids[.y]
          .y
        })
      } else if (!is.null(dim(.x))) {
        row.names(.x) <- ids
      } else if (length(.x) < n) {
        names(.x) <- ids[.x]
      } else {
        names(.x) <- ids
      }

      .x
    })
  }

  res$corr.mei <- corr.mei
  class(res) <- c("depthgram", class(res))
  res
}

#' @rdname depthgram
#' @export
depthgram.fData <- function(Data,
                            marginal_outliers = FALSE,
                            boxplot_factor = 1.5,
                            outliergram_factor = 1.5,
                            ids = NULL) {
  depthgram(
    Data = list(Data$values),
    marginal_outliers = marginal_outliers,
    boxplot_factor = boxplot_factor,
    outliergram_factor = outliergram_factor,
    ids = ids
  )
}

#' @rdname depthgram
#' @export
depthgram.mfData <- function(Data,
                             marginal_outliers = FALSE,
                             boxplot_factor = 1.5,
                             outliergram_factor = 1.5,
                             ids = NULL) {
  depthgram(
    Data = toListOfValues(Data),
    marginal_outliers = marginal_outliers,
    boxplot_factor = boxplot_factor,
    outliergram_factor = outliergram_factor,
    ids = ids
  )
}

#' Specialized method to plot 'depthgram' objects
#'
#' This function plots the three 'DepthGram' representations from the output of
#' the \code{\link{depthgram}} function.
#'
#' @param x An object of class `depthgram` as output by the
#'   \code{\link{depthgram}} function.
#' @param limits A boolean specifying whether the empirical limits for outlier
#'   detection should be drawn. Defaults to `FALSE`.
#' @param ids A character vector specifying labels for individual observations.
#'   Defaults to `NULL`, in which case observations will named by their id
#'   number in order of appearance.
#' @param print A boolean specifying whether the graphical output should be
#'   optimized for printed version. Defaults to `FALSE`.
#' @param plot_title A character string specifying the main title for the plot.
#'   Defaults to `""`, which means no title.
#' @param shorten A boolean specifying whether labels must be shorten to 15
#'   characters. Defaults to `TRUE`.
#' @param col Color palette used for the plot. Defaults to `NULL`, in which case
#'   a default palette produced by the \code{\link[grDevices]{hcl}} function is
#'   used.
#' @param pch Point shape. See \code{\link[plotly]{plotly}} for more details.
#'   Defaults to `19`.
#' @param sp Point size. See \code{\link[plotly]{plotly}} for more details.
#'   Defaults to `2`.
#' @param st Label size. See \code{\link[plotly]{plotly}} for more details.
#'   Defaults to `4`.
#' @param sa Axis title sizes. See \code{\link[plotly]{plotly}} for more
#'   details. Defaults to `10`.
#' @param text_labels A character vector specifying the labels for the
#'   individuals. It is overridden if `limits = TRUE`, for which only outliers
#'   labels are shown. See \code{\link[plotly]{plotly}} for more details.
#'   Defaults to `""`.
#' @param ... Other arguments to be passed to the base \code{\link[base]{plot}}
#'   function. Unused.
#'
#' @return A list with the following items:

#' - `p`: list with all the interactive (plotly) depthGram plots;
#' - `out`: outliers detected;
#' - `colors`: used colors for plotting.
#'
#' @references
#' Aleman-Gomez, Y., Arribas-Gil, A., Desco, M. Elias-Fernandez, A., and Romo,
#' J. (2021). "Depthgram: Visualizing Outliers in High Dimensional Functional
#' Data with application to Task fMRI data exploration".
#'
#' @export
#'
#' @examples
#' N <- 1e2
#' P <- 1e2
#' grid <- seq(0, 1, length.out = P)
#' Cov <- exp_cov_function(grid, alpha = 0.3, beta = 0.4)
#'
#' Data <- list()
#' Data[[1]] <- generate_gauss_fdata(
#'   N,
#'   centerline = sin(2 * pi * grid),
#'   Cov = Cov
#' )
#' Data[[2]] <- generate_gauss_fdata(
#'   N,
#'   centerline = sin(2 * pi * grid),
#'   Cov = Cov
#' )
#' names <- paste0("id_", 1:nrow(Data[[1]]))
#' DG <- depthgram(Data, marginal_outliers = TRUE, ids = names)
#' plot(DG)
plot.depthgram <- function(x,
                           limits = FALSE,
                           ids = NULL,
                           print = FALSE,
                           plot_title = "",
                           shorten = TRUE,
                           col = NULL,
                           pch = 19, sp = 2, st = 4, sa = 10,
                           text_labels = "",
                           ...) {
  n <- length(x$mei.mbd.d)

  type <- c(
    "Dimensions DepthGram",
    "Time DepthGram",
    "Time/Correlation DepthGram"
  )

  if (is.null(ids)) ids <- as.character(1:n)

  x <- data.frame(
    ID = rep(ids, 3),
    mei.mbd = c(x$mei.mbd.d, x$mei.mbd.t, x$mei.mbd.t2),
    mbd.mei = c(x$mbd.mei.d, x$mbd.mei.t, x$mbd.mei.t2),
    type = rep(type, each = n)
  )

  if (is.null(col)) {
    hues <- seq(15, 375, length = n + 1)
    color <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  } else
    color <- col

  out <- NULL

  if (limits) {
    P2 <- function(x, n) {
      a0 <- 2 / n
      a2 <- -n / (2 * (n - 1))
      a0 + x + a2*x^2
    }
    meis <- seq(0, 1, n)
    x$meis <- rep(meis, 3)
    x$par <- P2(x$meis, n)
    distp <- x$mbd.mei - P2(1 - x$mei.mbd, n)
    q3 <- stats::quantile(distp, 0.75)
    q1 <- stats::quantile(distp, 0.25)
    x$par2 <- x$par + q3 + 1.5 * (q3 - q1)
    out <- unique(which(distp > q3 + 1.5 * (q3 - q1)) %% n)
    pch <- rep(1, n) # empty circle
    pch[out] <- 19 # solid circle
    text_labels <- rep("", n)
    text_labels[out] <- ids[out]

    if (is.null(col)) {
      hues <- seq(15, 375, length = length(out) + 1)
      color.out <- grDevices::hcl(h = hues, l = 65, c = 100)[1:length(out)]
      color <- rep(8, n)
      color[out] <- color.out
    }
  }

  if (print) {
    sp <- 3
    st <- 5
    sa <- 14
  }

  plots <- list()

  for (i in 1:3) {
    dat <- x[which(x$type == type[i]), ]

    plots[[i]] <- ggplot(dat, aes(x = 1 - .data$mei.mbd, y = .data$mbd.mei)) +
      geom_point(aes(group = .data$ID), color = color, size = sp, shape = pch) +
      geom_text(
        label = text_labels,
        color = color,
        hjust = -0.15,
        vjust = -0.15,
        size = st
      ) +
      xlim(c(0, 1.005)) +
      ylim(c(0, 0.525)) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = sa),
        axis.text = element_text(size = sa - 2),
        title = element_text(size = sa)
      )

    if (limits) {
      plots[[i]] <- plots[[i]] +
        geom_line(aes(x = .data$meis, y = .data$par) , col = 1, na.rm = TRUE) +
        geom_line(aes(x = .data$meis, y = .data$par2), col = 1, na.rm = TRUE, lty = 2)
    }

    if (i == 1)
      pt <- plot_title
    else
      pt <- ""

    plots[[i]] <- plots[[i]] +
      facet_wrap(~ .data$type) +
      ggtitle(pt)
  }

  p <- plots %>%
    plotly::subplot(shareY = TRUE, shareX =TRUE) %>%
    plotly::layout(
      title = plot_title,
      yaxis = list(title = "MBD(MEI)"),
      xaxis = list(title = "1-MEI(MBD)")
    )

  print(p)

  list(
    p = list(
       dimDG = plots[[1]],
      timeDG = plots[[2]],
      corrDG = plots[[3]],
      fullDG = p
    ),
    out = out,
    color = color
  )
}
