#' @param sasPath single character string containing path or URL for some SAS
#'   code. The source SAS code is stored in extdata or at the location specified
#'   in URL. In development mode in \code{devtools}, this is actually
#'   inst/extdata. This function is internal, since it is used to generate data
#'   which ends up in the distributed package. However, the package user can
#'   verify that the code creates the distributed R data. Also, changes to the
#'   original SAS code can be used to regenerate the R data by the user, without
#'   waiting for a package release.
