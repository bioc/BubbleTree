
#' @title saveXLS
#' @rdname saveXLS
#' @aliases saveXLS
#' @param dat.lst dataframe
#' @param xls.fn filename
#' @param row.names row names
#' @param ... misc
#' @return new Excel file
#' @example examples/saveXLS-Ex.R
saveXLS <- function(dat.lst, xls.fn, row.names=FALSE, ...) {
    warning("Be cautious of embedded quotes 
            (using remove.embedded.quotes)!!!\n")

    if(row.names) {
        dat.lst <- lapply(dat.lst, function(x) cbind(row.names(x), x))
    }

    with(dat.lst, {
        WriteXLS(names(dat.lst), xls.fn)
    })
}

#' @title btcompare
#' @rdname btcompare
#' @aliases btcompare
#' @param set1 first set
#' @param set2 second set to compare
#' @return combined, unique list of genes
#' @example examples/btcompare-Ex.R
btcompare <- function(set1, set2) {
    comm <- intersect(set1, set2)
    uniq.1 <- setdiff(set1, comm)
    uniq.2 <- setdiff(set2, comm)
    cat(sprintf("%d are common, whereas %d and %d are unique in each dataset\n",
                length(comm),
                length(uniq.1),
                length(uniq.2)))

    return(list(comm=comm, uniq.1=uniq.1, uniq.2=uniq.2))
}
