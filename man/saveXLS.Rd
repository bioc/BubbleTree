% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{saveXLS}
\alias{saveXLS}
\title{saveXLS}
\usage{
saveXLS(dat.lst, xls.fn, row.names = FALSE, ...)
}
\arguments{
\item{dat.lst}{dataframe}

\item{xls.fn}{filename}

\item{row.names}{row names}

\item{...}{misc}
}
\value{
new Excel file
}
\description{
saveXLS
}
\examples{

load(system.file("data", "allCall.lst.RData", package="BubbleTree"))

all.summary <- plyr::ldply(allCall.lst, function(.Object) {
    purity <- .Object@result$prev[1]
    adj <- .Object@result$ploidy.adj["adj"]
    # when purity is low the calculation result is not reliable
    ploidy <- (2*adj -2)/purity + 2  
    
    with(.Object@result,
         return(c(Purity=round(purity,3),
                  Prevalences=paste(round(prev,3), collapse=", "),
                  "Tumor ploidy"=round(ploidy,1))))
}) \%>\% plyr::rename(c(".id"="Sample"))

xls.filename <- paste("all_summary", "xlsx", sep=".")
saveXLS(list(Summary=all.summary), xls.filename)
}
