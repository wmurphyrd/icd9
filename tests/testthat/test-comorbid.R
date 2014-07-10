context("comorbidities, slow tests")

#randomSampleAhrq <- sample(unname(c(ahrqComorbid, recursive = TRUE)), replace = TRUE, size = n)
testSimplePatients <- data.frame(
  visitId = as.character(c(1000, 1000, 1000, 1001, 1001, 1002)),
  icd9 = c("27801", "7208", "25001", "34400", "4011", "4011"),
  poa = factor(c("Y", "N", "Y", "N", "Y", "N")),
  stringsAsFactors = FALSE
)

n = 100
np = 2
set.seed(1441)
testRandomPatients <- data.frame(
  visitId = as.character(sample(seq(1, np), replace = TRUE, size = n)),
  #icd9 = as.character(floor(runif(min = 1000, max = 99999, n = n))),
  icd9 = as.character(sample(unname(unlist(ahrqComorbid)), size = n)),
  poa = as.factor(sample(x=c("Y", "N", "n", "n", "y", "X", "E", "", NA), replace = TRUE, size = n)),
  stringsAsFactors = FALSE
)

# not real data, but designed to be similar
testRealPatients <- structure(
  list(visitId = as.character(c(207210584L, 207210584L, 207210584L,
                                207210584L, 207210584L, 207210600L, 207210600L,
                                207210600L, 207210600L, 207210600L, 207210600L,
                                207210600L, 207210600L, 207210600L, 207210600L,
                                207210600L, 207210600L, 207210600L, 207210618L, 207210618L)),
       icd9 = c("19907", "V4346", "26189", "19318", "20149", "30434", "14291",
                "28012", "24940", "3431", "20971", "20040", "24948", "28063",
                "16876", "30048", "16502", "1797", "16016", "20914"),
       poa = c("N", "N", "N", "Y", "Y", "Y", "Y", "Y", "Y", "Y",
               "Y", "Y", "Y", "Y", "E", "E", "Y", "Y", "Y", "N")),
  .Names = c("visitId", "icd9", "poa"),
  class = "data.frame"
)

#TODO: tests for these:
#testMultiDzOnePat
#testMultiDzMultiPats
#testEmptyInvalidPats

#testRagged


# test that the C++ version always give the same results as the R version:
# doesn't mean they're right, of course.
for (p in c("testSimplePatients", "testRandomPatients", "testRealPatients")) {
  pts <- get(p)
  test_that(paste("R and C++ gives same result with ", p), {
    expect_identical(
      icd9Comorbidities(visitId = pts$visitId, icd9 = pts$icd9, icd9Mapping = ahrqComorbid, doWithCpp = TRUE),
      icd9Comorbidities(visitId = pts$visitId, icd9 = pts$icd9, icd9Mapping = ahrqComorbid, doWithCpp = FALSE),
      info = p)
  })
}

# test the c++ and R versions against each other
test_that(paste("icd9 comorbidities extra id/code to mismatch length"), {
  expect_error(icd9Comorbidities(pts$visitId, c(pts$icd9, "12345"), ahrqComorbid, doWithCpp = TRUE))
  expect_error(icd9Comorbidities(c(pts$visitId, "123"), pts$icd9, ahrqComorbid, doWithCpp = TRUE))
  expect_error(icd9Comorbidities(pts$visitId, c(pts$icd9, "12345"), ahrqComorbid, doWithCpp = FALSE))
  expect_error(icd9Comorbidities(c(pts$visitId, "123"), pts$icd9, ahrqComorbid, doWithCpp = FALSE))
})

test_that(paste("mapping is not a list with"), {
  expect_error(icd9Comorbidities(pts$visitId, pts$icd9, "ceci n'est pas un mapping", doWithCpp = FALSE))
  expect_error(icd9Comorbidities(pts$visitId, pts$icd9, "ceci n'est pas un mapping", doWithCpp = TRUE))
})

test_that(paste("data frame not allowed anymore with"), {
  expect_error(icd9Comorbidities(testRealPatients, ahrqComorbid, doWithCpp = FALSE))
  expect_error(icd9Comorbidities(testRealPatients, ahrqComorbid, doWithCpp = TRUE))
  expect_error(icd9Comorbidities(testRealPatients, testRealPatients, ahrqComorbid, doWithCpp = FALSE))
  expect_error(icd9Comorbidities(testRealPatients, testRealPatients, ahrqComorbid, doWithCpp = TRUE))
})

test_that(paste("non-character icd9 input with"), {
  # TODO: these get casted to charactervector silently!
  #expect_error(icd9Comorbidities(c("1", "2", "3"), c(1,2,3))))
  #expect_error(icd9Comorbidities(c("1", "2", "3"), c(1L,2L,3L))))
  # todo: what to do with factors?
})

test_that(paste("icd9 comorbidities using"), {

  for (pts in c(testSimplePatients, testRandomPatients, testRealPatients)) {
    p <- get(pts)
    for (tf in c(TRUE, FALSE)) {
      m <- icd9Comorbidities(visitId = p$visitId, icd9 = p$icd9, icd9Mapping = ahrqComorbid, doWithCpp = tf)
      expect_equal(colnames(m), names(ahrqComorbid), info = p)
      expect_equal(sort(unique(asCharacterNoWarn(p$visitId))), sort(rownames(m)), info = p) # should match exactly, but no order guarantee?
      expect_is(m, "matrix", info = p)
      expect_equal(typeof(m), "logical", info = p)
    } # end for T/F
  } # end for test datasets
})


# ptdflogical <- logicalToBinary(ptdf)  # TODO: test elsewhere
#expect_true(all(sapply(names(ahrqComorbid), function(x) class(ptdflogical[, x])) == "integer"))

#   expect_equal(
#     logicalToBinary(data.frame(a = c("jack", "hayley"), b = c(TRUE, FALSE), f = c(TRUE, TRUE))),
#     data.frame(a = c("jack", "hayley"), b = c(1, 0), f = c(1, 1))
#   )

test_that("no NA values in the co-morbidity lists", {
  expect_false(any(is.na(unlist(unname(ahrqComorbid)))))
  expect_false(any(is.na(unlist(unname(ahrqComorbidAll)))))
  expect_false(any(is.na(unlist(unname(quanDeyoComorbid)))))
  expect_false(any(is.na(unlist(unname(quanElixhauserComorbid)))))
  expect_false(any(is.na(unlist(unname(elixhauserComorbid)))))
})

test_that("can convert mappings to and from short and decimal", {
  # TODO: decaimal to short and back
  expect_equal(icd9MappingDecimalToShort(icd9MappingShortToDecimal(ahrqComorbid)), ahrqComorbid)
})

test_that("built-in icd9 to comorbidity mappings are all valid", {
  expect_true(icd9ValidMappingShort(ahrqComorbid))
  expect_true(icd9ValidMappingShort(quanDeyoComorbid))
  expect_true(icd9ValidMappingShort(quanElixhauserComorbid))
  expect_true(icd9ValidMappingShort(elixhauserComorbid))
})

test_that("ahrq icd9 mappings are all generated from the current generation code", {
  # same but from source data. Should be absolutely identical.
  expect_identical(ahrqComorbid, parseAhrqSas(condense = FALSE, save = FALSE, returnAll = FALSE))
  # same but from source data. Should be absolutely identical.
  expect_identical(ahrqComorbidAll, parseAhrqSas(condense = FALSE, save = FALSE, returnAll = TRUE))
})

test_that("quan charlson icd9 mappings are all generated from the current generation code", {
  expect_identical(quanDeyoComorbid, parseQuanDeyoSas(condense = FALSE, save = FALSE))
})

test_that("quan elixhauser icd9 mappings are all generated from the current generation code", {
  expect_identical(quanElixhauserComorbid, parseQuanElixhauser(condense = FALSE, save = FALSE))
})

test_that("elixhauser icd9 mappings are all generated from the current generation code", {
  expect_identical(elixhauserComorbid, parseElixhauser(condense = FALSE, save = FALSE))
})

test_that("can condense the big lists of comorbidities without errors", {
  # this is a useful test because the data weren't generated by just expanding
  # base ranges (which is how the condense works in reverse)
  ahrq <- lapply(ahrqComorbid, icd9CondenseShort)
  quanDeyo <- lapply(quanDeyoComorbid, icd9CondenseShort)
  quanElixhauser <- lapply(quanElixhauserComorbid, icd9CondenseShort)
  elixhauser <- lapply(elixhauserComorbid, icd9CondenseShort)
  expect_is(ahrq, class = "list")
  expect_is(elixhauser, class = "list")
  expect_is(quanDeyo, class = "list")
  expect_is(quanElixhauser, class = "list")
  # the comorbidity mappings save in \code{data} should not be condensed.
  expect_false(identical(ahrq, ahrqComorbid))
  expect_false(identical(elixhauser, elixhauserComorbid))
  expect_false(identical(quanDeyo, quanDeyoComorbid))
  expect_false(identical(quanElixhauser, quanElixhauserComorbid))
})

test_that("ahrq make sure all the children are listed in the saved data.", {
  ahrq <- lapply(ahrqComorbid, icd9ChildrenShort)
  expect_equal(ahrq, ahrqComorbid)
})
test_that("elixhauser make sure all the children are listed in the saved data.", {
  elixhauser <- lapply(elixhauserComorbid, icd9ChildrenShort)
  expect_equal(elixhauser, elixhauserComorbid)
})
test_that("quan charlson make sure all the children are listed in the saved data.", {
  quanDeyo <- lapply(quanDeyoComorbid, icd9ChildrenShort)
  expect_equal(quanDeyo, quanDeyoComorbid)
})
test_that("quan elixhauser make sure all the children are listed in the saved data.", {
  quanElixhauser <- lapply(quanElixhauserComorbid, icd9ChildrenShort)
  expect_equal(quanElixhauser, quanElixhauserComorbid)
})

test_that("condense an ICD-9 code set to minimal group", {
  expect_equal(sort(icd9CondenseShort("98799" %i9s% "98901")), sort(c("98799", "988", "98900", "98901")))
  # TODO: more tests
})


# test_that("AHRQ interpretation at least returns something reasonable", {
#   result <- parseAhrqSas(sasPath = system.file("extdata", "comformat2012-2013.txt", package="icd9"), save = FALSE)
#   expect_that(result, is_a("list"))
#   expect_true(length(result) > 10)
# })

test_that("HTN subgroups all worked", {
  # pick one subcategory
  expect_true(all(ahrqComorbidAll$HTNPREG %in% ahrqComorbid$HTNCX))

  # and we didn't drop any:
  expect_true(all(ahrqComorbidAll$HTNCX %in% ahrqComorbid$HTNCX))
  expect_true(all(ahrqComorbidAll$CHF %in% ahrqComorbid$CHF))
  expect_true(all(ahrqComorbidAll$RENLFAIL %in% ahrqComorbid$RENLFAIL))

})
