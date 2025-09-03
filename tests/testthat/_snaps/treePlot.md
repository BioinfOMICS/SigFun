# treePlot handles errors correctly

    Code
      treePlot(GSE181574.sigfun, showCategory = 2)
    Condition
      Error in `treePlot()`:
      ! `showCategory` must large than "5".

---

    Code
      treePlot(GSE181574.sigfun, showCategory = c("GOMF_STEROID_HYDROXYLASE_ACTIVITY",
        "WP_FATTY_ACID_OMEGAOXIDATION"))
    Condition
      Error in `treePlot()`:
      ! `showCategory` must provide at least "5" pathway names.

