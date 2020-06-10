## Test environments
I used the Github Action script
[check-standard](https://github.com/r-lib/actions/blob/master/examples/check-standard.yaml)
provided by the RStudio team to setup automatic environment checks. This
workflow runs `R CMD check` via the **rcmdcheck** package on the three major OSs
(linux, macOS and Windows) with the current release version of R, and R-devel.

## R CMD check results
There were no ERRORs.

There was 1 WARNING when running `R CMD check` locally because I have not
installed **qpdf** on my machine:

```
* checking data for ASCII and uncompressed saves ... OK
   WARNING
  ‘qpdf’ is needed for checks on size reduction of PDFs
```

There were no NOTEs.
