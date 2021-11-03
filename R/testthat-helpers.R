# To use expect_snapshot_file() you'll typically need to start by writing
# a helper function that creates a file from your code, returning a path
save_png <- function(code, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  code

  path
}

# You'd then also provide a helper that skips tests where you can't
# be sure of producing exactly the same output
expect_snapshot_plot <- function(name, code) {
  # Other packages might affect results
  # skip_if_not_installed("ggplot2", "2.0.0")
  # Or maybe the output is different on some operation systems
  # skip_on_os("windows")
  # You'll need to carefully think about and experiment with these skips
  # skip_on_ci()

  name <- paste0(name, ".png")

  # Announce the file before touching `code`. This way, if `code`
  # unexpectedly fails or skips, testthat will not auto-delete the
  # corresponding snapshot file.
  testthat::announce_snapshot_file(name = name)

  path <- save_png(code)
  testthat::expect_snapshot_file(path, name)
}
