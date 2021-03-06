% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports.R, R/comorbid.R
\name{icd9InReferenceCode}
\alias{\%i9in\%}
\alias{icd9InReferenceCode}
\title{match ICD9 codes}
\usage{
icd9InReferenceCode(icd9, icd9Reference, isShort, isShortReference = TRUE)

icd9 \%i9in\% icd9Reference
}
\arguments{
\item{isShort}{single logical value which determines whether the ICD-9 code
provided is in short (TRUE) or decimal (FALSE) form. Where reasonable, this
is guessed from the input data.}

\item{isShortReference}{logical, see argument \code{isShort}}

\item{"icd9,icd9Reference"}{is
a character vector or factor of ICD-9 codes. If fewer than five characters
is given in a code, then the digits are greedily assigned to hundreds, then
tens, then units, before the decimal parts. E.g. "10" becomes "010", not
"0010".}
}
\value{
logical vector
}
\description{
Finds children of \code{icd9Reference} and looks for \code{icd9} in the
  resulting vector.
}
\examples{
#\%i9in\% assumes both test code(s) and reference set of codes are \\emph{short}
"1024" \%i9in\% "102"
"1024" \%i9in\% c("102","1025")
c("102", "1024","1025") \%i9in\% "102"
c("102", "1024","1025") \%i9in\% c("1024", "1025")
c("102", "1024","1025") \%i9in\% c("102", "1024", "1025")
}
\keyword{internal}

