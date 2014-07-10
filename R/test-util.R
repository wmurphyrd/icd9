# internal test-data generation functions may be better to leave these in the
# tests folder, but would like the convenience of calling them from the internal
# package namespace when in dev_mode

randomPatients <- function(n = 50000, np = 20) {
  pts <- round(n/np)
  data.frame(
    visitId = sample(seq(1, pts), replace = TRUE, size = n),
    icd9 = randomShortIcd9(n),
    poa = as.factor(
      sample(x = c("Y","N", "n", "n", "y", "X","E","",NA), replace = TRUE, size = n))
  )
}

randomShortIcd9 <- function(n = 50000)
  as.character(floor(runif(min=1, max = 99999, n=n))) # tolerate <3 digits?

randomDecimalIcd9 <- function(n = 50000)
  paste(
    round(runif(min = 1, max = 999, n = n)),
    sample(icd9ExpandMinor(), replace = TRUE, size = n),
    sep = "."
  )
