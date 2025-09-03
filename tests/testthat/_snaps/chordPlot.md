# chordPlot handles fontSize errors correctly

    Code
      chordPlot(GSE181574.sigfun, fontSize = "invalidSize")
    Condition
      Error in `chordPlot()`:
      ! The input `fontSize` must be a <numeric>.

# chordPlot handles fontSize boundary errors correctly

    Code
      chordPlot(GSE181574.sigfun, fontSize = 0)
    Condition
      Error in `chordPlot()`:
      ! The input `fontSize` must be a numeric value between 0 and 1.

---

    Code
      chordPlot(GSE181574.sigfun, fontSize = 1.5)
    Condition
      Error in `chordPlot()`:
      ! The input `fontSize` must be a numeric value between 0 and 1.

# chordPlot handles invalid pathway errors correctly

    Code
      chordPlot(GSE181574.sigfun, showCategory = c("invalid_pathway1",
        "invalid_pathway2"))
    Condition
      Error in `.extractGeneSets()`:
      ! The value provided to `showCategory` cannot be found in SE_data.fgsea@result.

