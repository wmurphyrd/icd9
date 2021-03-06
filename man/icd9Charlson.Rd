% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/score.R
\name{icd9Charlson}
\alias{icd9Charlson}
\alias{icd9Charlson.data.frame}
\alias{icd9CharlsonComorbid}
\title{Calculate Charlson Comorbidity Index (Charlson Score)}
\usage{
icd9Charlson(x, visitId = NULL, return.df = FALSE,
  stringsAsFactors = getOption("stringsAsFactors"), ...)

\method{icd9Charlson}{data.frame}(x, visitId = NULL, return.df = FALSE,
  stringsAsFactors = getOption("stringsAsFactors"), ...)

icd9CharlsonComorbid(x, visitId = NULL, applyHierarchy = FALSE)
}
\arguments{
\item{x}{data frame containing a column of visit or patient identifiers, and
a column of ICD-9 codes. It may have other columns which will be ignored.
By default, the first column is the patient identifier and is not counted.
If \code{visitId} is not specified, the first column is used.}

\item{visitId}{The name of the column in the data frame which contains the
patient or visit identifier. Typically this is the visit identifier, since
patients come leave and enter hospital with different ICD-9 codes. It is a
character vector of length one. If left empty, or \code{NULL}, then an
attempt is made to guess which field has the ID for the patient encounter
(not a patient ID, although this can of course be specified directly). The
guesses proceed until a single match is made. Data frames may be wide with
many matching fields, so to avoid false positives, anything but a single
match is rejected. If there are no successful guesses, and \code{visitId}
was not specified, then the first column of the data frame is used.}

\item{return.df}{single logical value, if true, a two column data frame will
be returned, with the first column named as in input data frame (i.e.
\code{visitId}), containing all the visits, and the second column
containing the Charlson Comorbidity Index.}

\item{stringsAsFactors}{single logical, passed on when constructing
data.frame if \code{return.df} is \code{TRUE}. If the input data frame
\code{x} has a factor for the visitId, this is not changed, but a
non-factor visitId may be converted or not converted according to your
system default or this setting.}

\item{...}{further arguments to pass on to \code{icd9ComorbidQuanDeyo}, e.g.
\code{icd9Field}}

\item{applyHierarchy}{single logical value, default is FALSE. If TRUE, will
drop DM if DMcx is present, etc.}
}
\description{
Charlson score is calculated in the basis of the Quan revision
  of Deyo's ICD-9 mapping. (Peptic Ulcer disease no longer warrants a point.)
  Quan published an updated set of scores, but it seems most people use the
  original scores for easier comaprison between studies, even though Quan's
  were more predictive.
}
\details{
Per Quan, "The following comorbid conditions were mutually
  exclusive: diabetes with chronic complications and diabetes without chronic
  complications; mild liver disease and moderate or severe liver disease; and
  any malignancy and metastatic solid tumor.""
}
\section{Methods (by class)}{
\itemize{
\item \code{data.frame}: Charlson scores from data frame of visits and ICD-9 codes
}}
\examples{
mydf <- data.frame(visitId = c("a", "b", "c"),
                   icd9 = c("441", "412.93", "044.9"))
cmb <- icd9ComorbidQuanDeyo(mydf, isShort = FALSE, applyHierarchy = TRUE)
cmb
icd9Charlson(mydf, isShort = FALSE)
icd9Charlson(mydf, isShort = FALSE, return.df = TRUE)
icd9CharlsonComorbid(cmb)
}

