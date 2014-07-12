#' benchmark and profile major functions with larger data sets
#'
#' \code{icd9} package is intended to be used with large data sets, with
#' millions or rows. Performance of core functions is therefore of some
#' importance, after correctness. R package test code is for correctness,
#' whereas this script stresses the core functions, and looks for bottlenecks.

icd9Benchmark <- function() {
  require(microbenchmark)
  # generate large data set: this is copied from test-ICD9.R for now...
  set.seed(1441)
  n <- 1E7 # 10 million rows
  np <- 20 # 20 icd9 codes per patient

  rpts <- randomPatients(n)

  tmp <- tempfile(fileext = ".Rprof")
  Rprof(filename = tmp, line.profiling = T, memory.profiling = T)
  capture.output(icd9Comorbidities(rpts))
  Rprof(NULL)

  #summaryRprof(filename = tmp, memory = "stats", lines = "both")
  summaryRprof(filename = tmp, memory = "both", lines = "show")

  #microbenchmark(times = 10, icd9ShortToParts(randomShortIcd9(5E+5)))
  #microbenchmark(times = 10, icd9ShortToPartsSlow(randomShortIcd9(5E+5)))
  #microbenchmark(times = 10, icd9ShortToPartsList(randomShortIcd9(5E+5)))

  microbenchmark(times = 50, trim(randomShortIcd9))
  microbenchmark(times = 50, strip(randomShortIcd9))

  # initializing empty data frame
  microbenchmark(data.frame(matrix(ncol = 2, nrow = 100000)))
  microbenchmark(data.frame(major = character(100000), minor = character(100000)))

}

benchmarkComorbid <- function(n = 1e4, np = 5, times = 10) {
  require(microbenchmark)
  patientData <- randomPatients(n = n, np = np)
  print(
    microbenchmark(
      icd9ComorbiditiesLongR(patientData$visitId, patientData$icd9, icd9Mapping = ahrqComorbid),
      icd9ComorbiditiesLongCpp(patientData$visitId, patientData$icd9, icd9Mapping = ahrqComorbid),
      times = times
    )
  )
}

# the overhead is all in the bloated \code{icd9InReferenceCode} function. This
# is ameliorated by memoisation, at the cost of memory usage and complexity.
# over two million rows, the three methods (with memoisation) are the same.
# Unit: seconds
# expr      min       lq   median       uq      max neval
# one 13.26761 13.70637 13.75466 14.09683 14.28764    10
# two 13.70615 13.77829 13.85725 13.96494 14.02498    10
# three 13.21530 13.47809 13.78816 13.99146 14.58961    10

benchmarkComorbidRmethod <- function() {
  require("microbenchmark")

  tm = 1e6
  icd9df <- data.frame(
    visitId = rep(c("patcom100", "patcom222"), times = tm),
    icd9Short = rep(c("4280", "3979"), times = tm),
    poa = rep(c("Y", "N"), times = tm),
    stringsAsFactors = FALSE
  )
  icd9Mapping <- ahrqComorbid

  # someday/maybe: do with memoisation:
  #require(memoise)
  #memSpawnReferenceChildren <- memoise::memoise(icd9SpawnReferenceChildren)

  microbenchmark(times = 10L,
    one = {
    cbind(
      visitId = icd9df[['visitId']],
      matrix(
        dimnames = list(icd9df[['visitId']], names(icd9Mapping)), # force matrix so that n=1 doesn't give stupid output.
        nrow = length(icd9df[['icd9Short']]),
        data = vapply(
          X = names(icd9Mapping),
          FUN.VALUE = rep(FALSE, length(icd9df[['icd9Short']])),
          FUN = function(comorbidity) {
            icd9InReferenceCode(
              # drop factor down to character codes #TODO: is this necessary or desirable?
              asCharacterNoWarn(icd9df[['icd9Short']]),
              # provide vector of base ICD9 codes for this comorbidity group
              icd9Mapping[[comorbidity]]
            )
          }
        )
      )
    )
  },
  two = {
    ic <- asCharacterNoWarn(icd9df[['icd9Short']])
    cbind(
      visitId = icd9df[['visitId']],
      matrix(
        dimnames = list(icd9df[['visitId']], names(icd9Mapping)), # force matrix so that n=1 doesn't give stupid output.
        nrow = length(icd9df[['icd9Short']]),
        data = vapply(
          X = names(icd9Mapping),
          FUN.VALUE = rep(FALSE, length(icd9df[['icd9Short']])),
          FUN = function(comorbidity) {
            icd9InReferenceCode(
              # drop factor down to character codes #TODO: is this necessary or desirable?
              ic,
              # provide vector of base ICD9 codes for this comorbidity group
              icd9Mapping[[comorbidity]]
            )
          }
        )
      )
    )
  },
  three = {
    ic <- asCharacterNoWarn(icd9df[['icd9Short']])
    cbind(
      visitId = icd9df[['visitId']],
      matrix(
        dimnames = list(icd9df[['visitId']], names(icd9Mapping)), # force matrix so that n=1 doesn't give stupid output.
        nrow = length(icd9df[['icd9Short']]),
        data = vapply(
          X = names(icd9Mapping),
          FUN.VALUE = rep(FALSE, length(icd9df[['icd9Short']])),
          FUN = function(comorbidity) ic %in% icd9Mapping[[comorbidity]]
        )
      )
    )
  }
  )
}
