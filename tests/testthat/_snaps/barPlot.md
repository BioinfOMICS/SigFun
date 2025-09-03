# barPlot handles errors correctly

    Code
      barPlot(GSE181574.sigfun, color = "invalid_method")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pvalue", "p.adjust", "NES"

---

    Code
      barPlot(GSE181574.sigfun, showCategory = c("invalid_pathway1",
        "invalid_pathway2"))
    Condition
      Error in `.extractGeneSets()`:
      ! The value provided to `showCategory` cannot be found in SE_data.fgsea@result.

