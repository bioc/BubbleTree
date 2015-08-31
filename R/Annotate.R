#' @title Annotate
#' @description Annotate
#' @docType package
#' @name Annotate
#' @examples 
#' annotate <- new("Annotate")
utils::globalVariables(c("."))

gr <- setClass("GenomicRanges")

Annotate <- setClass(
    "Annotate",

    representation(chr="character",
                   beg="numeric",
                   end="numeric",
                   critical.genes="character"),

    prototype = list(chr="",
                     beg=1,
                     end=500,
                     critical.genes=NULL)
)

setMethod("initialize",
          "Annotate",
          function(.Object, chr="", beg=0, end=0, critical.genes="") {
              .Object@chr <- character()
              .Object@beg <- numeric()
              .Object@end <- numeric()
              .Object@critical.genes <- character()
              #.Object@gene.uni.clean.gr <- numeric()
              #.Object@cyto.gr <- numeric()
              .Object
          }
)

#' @export
#' @docType methods
#' @rdname annoByGenesAndCyto
setGeneric(name="annoByGenesAndCyto",
           def=function(.Object, chr, beg, end, 
                        critical.genes, gene.uni.clean.gr, cyto.gr) {
               standardGeneric("annoByGenesAndCyto")
           }
)

#' @title annoByGenesAndCyto
#' @rdname annoByGenesAndCyto
#' @aliases annoByGenesAndCyto
#' @param .Object the objet
#' @param chr the chromosome
#' @param beg genomic start coord
#' @param end genomic end coord
#' @param critical.genes set of critical genes
#' @param gene.uni.clean.gr gr object of genes
#' @param cyto.gr gr object of cyto positions
#' @return list of annotation for genes and cytobands
#' @example examples/annoByGenesAndCyto-Ex.R
setMethod("annoByGenesAndCyto",
          "Annotate",
          function(.Object, chr, beg, end, 
                   critical.genes, gene.uni.clean.gr, cyto.gr) {
              ann.grs <- gene.uni.clean.gr
              if(!is.null(critical.genes)) {
                  ann.grs <- 
                      gene.uni.clean.gr[
                          gene.uni.clean.gr$gene.symbol %in% critical.genes ]
              }
              ann <- annoByOverlap(.Object, chr, beg, end, gene.grs=ann.grs)
              rv <- list()
              rv$ann <-  sapply(ann, paste, collapse=", ")

              cyto <- annoByOverlap(.Object, 
                                    chr, 
                                    beg, 
                                    end, 
                                    cyto.gr, 
                                    ann=cyto.gr$name)
              
              rv$cyto <- paste(sub("chr", "", chr), 
                               sapply(cyto, function(x) 
                                   paste(x[1], tail(x, 1), sep='-')), 
                               sep="")
              return(rv)
          }
)

setGeneric(name="annoByOverlap",
           def=function(.Object, chr, beg, end, 
                        gene.grs=NULL, gr=NULL, ann=gene.grs$gene.symbol) {
               standardGeneric("annoByOverlap")
           }
)

setMethod("annoByOverlap",
          "Annotate",
          function(.Object, chr, beg, end, 
                   gene.grs=NULL, gr=NULL, ann=gene.grs$gene.symbol) {

              if(is.null(gr)){
                  gr <- GRanges(chr, IRanges(beg, end))
              }
              ol <- as.data.frame(findOverlaps(gr, gene.grs))
              ol$ann <- ann[ol$subjectHits]
              rv <- dlply(ol, .(queryHits), function(df) df$ann)
              out <- list()
              out[1:length(gr)] <- ""
              out[as.numeric(names(rv))] <- rv
              return(out)
          }
)
