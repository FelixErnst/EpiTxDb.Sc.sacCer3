#' @title Annotation package for \code{EpiTxDb} objects
#'
#' @author Felix G M Ernst [aut]
#'
#' @description
#' This package loads one or more \code{EpiTxDb} objects. Such \code{EpiTxDb}
#' objects are an R interface to prefabricated databases contained by this
#' package.
#' 
#' The names of any objects exposed by this package indicate the origin and
#' resources exposed. So for example \code{EpiTxDb.Sc.sacCer3.tRNAdb} would be a
#' \code{EpiTxDb} object for Saccharomyces cerevisia data from tRNAdb build
#' based on the sacCer3 genome build.
#' 
#' @return a \code{\link[EpiTxDb:EpiTxDb-class]{EpiTxDb}} object 
#' 
#' @seealso
#' \itemize{
#' \item{\code{\link[EpiTxDb:modifications]{modifications}}}
#' \item{\code{\link[EpiTxDb:modifications]{reactions}}}
#' \item{\code{\link[EpiTxDb:modifications]{specifies}}}
#' \item{\code{\link[EpiTxDb:modifications]{modificationsByTranscript}}}
#' \item{\code{\link[EpiTxDb:modifications]{modifiedSeqsByTranscript}}}
#' }
#' 
#' @docType package
#' @name EpiTxDb.Sc.sacCer3
#' 
#' @examples 
#' EpiTxDb.Sc.sacCer3.RMBase
NULL

#' @import EpiTxDb
NULL

.onLoad <- function(libname, pkgname)
{
    ns <- asNamespace(pkgname)
    path <- system.file("extdata", package = pkgname, lib.loc = libname)
    files <- dir(path)
    files <- files[grepl("*\\.sqlite",files)]
    for(i in seq_len(length(files))){
        db <- 
          AnnotationDbi::loadDb(system.file("extdata", files[[i]],
                                            package = pkgname,
                                            lib.loc = libname),
                                    packageName = pkgname)
        objname <- sub(".sqlite$", "", files[[i]])
        assign(objname, db, envir = ns)
        namespaceExport(ns, objname)
    }
}


#' @rdname EpiTxDb.Sc.sacCer3
#' @keywords datasets
#' 
#' @usage EpiTxDb.Sc.sacCer3.tRNAdb
"EpiTxDb.Sc.sacCer3.tRNAdb"
#' @rdname EpiTxDb.Sc.sacCer3
#' @usage EpiTxDb.Sc.sacCer3.RMBase
"EpiTxDb.Sc.sacCer3.RMBase"
