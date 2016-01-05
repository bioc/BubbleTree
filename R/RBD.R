#' @title RBD
#' @description RBD
#' @docType package
#' @name RBD
#' @examples 
#' rbd <- new("RBD")
utils::globalVariables(c("num.mark", "freq", "cnv.start", "cnv.end", 
                         "seg.mean", "."))

setClass("GenomicRanges")

RBD <- setClass(
    "RBD",

    representation(unimodal.kurtosis="numeric"),
    
    prototype(unimodal.kurtosis=-0.1)

)

setMethod("initialize",
          "RBD",
          function(.Object, unimodal.kurtosis=-0.1, ...) {
              .Object@unimodal.kurtosis <- unimodal.kurtosis
              .Object
          }
)

#' @export
#' @docType methods
#' @rdname makeRBD
setGeneric(name="makeRBD",
           def=function(.Object, unimodal.kurtosis=-0.1, ...) {
               standardGeneric("makeRBD")
           }
)

#' @title makeRBD
#' @description make the RBD object
#' @rdname makeRBD
#' @aliases makeRBD
#' @param .Object the object
#' @param unimodal.kurtosis kurtosis
#' @param snp.gr SNP GenomicRanges object
#' @param cnv.gr CNV GenomicRanges object
#' @param ... other input (not needed)
#' @return RBD object
#' @example examples/makeRBD-Ex.R
setMethod("makeRBD",
          "RBD",
          function(.Object, snp.gr, cnv.gr, unimodal.kurtosis=-0.1) {

              cs <- mergeSnpCnv(.Object, snp.gr, cnv.gr)
              
              cs1 <- cs %>% 
                  filter(!is.na(num.mark) & !is.na(seg.id) & !is.na(freq)) %>% 
                  group_by(seg.id,
                           seqnames,
                           cnv.start,
                           cnv.end,
                           num.mark,
                           seg.mean) %>% 
                  dplyr::summarise(kurtosis=e1071::kurtosis(freq, type=1),
                                   hds = ifelse(kurtosis >= unimodal.kurtosis,
                                                abs(median(freq, 
                                                           na.rm=TRUE)-0.5),
                                                median(abs(freq-0.5), 
                                                       na.rm=TRUE)),
                                   hds.sd= sd(abs(freq-0.5), na.rm=TRUE),
                                   het.cnt = length(freq) ) %>% 
                  plyr::rename(c("cnv.start"="start",
                                 cnv.end="end",
                                 seg.mean="lrr")) %>% as.data.frame

              # add those segment with low lrr < -1.25
              ids <- setdiff((1:length(cnv.gr))[cnv.gr$seg.mean < -1.25], 
                             cs1$seg.id)

              if(length(ids) >0) {
                  cs2 <- cnv.gr %>% 
                      as.data.frame %>% 
                      mutate(seg.id=1:length(cnv.gr)) %>% 
                      filter(seg.id %in% ids) %>% 
                      mutate(kurtosis = NA,
                             lrr=seg.mean,
                             hds=NA,
                             hds.sd=NA,
                             het.cnt=1)  %>% `[`(names(cs1))
                  cs1 <- rbind(cs1, cs2)
              }
              
              rbd <- with(cs1, GRanges(seqnames, IRanges(start, end)))
              elementMetadata(rbd) <- cs1[, ! names(cs1) %in% c("seqnames",
                                                                "start",
                                                                "end")]

              # calculate seg.size
              total.mark <- sum(rbd$num.mark, na.rm=TRUE)
              rbd$seg.size <- rbd$num.mark / total.mark * 100

              return(as.data.frame(rbd))
          }
)

#' @docType methods
#' @rdname mergeSnpCnv
setGeneric(name="mergeSnpCnv",
           def=function(.Object, snp.gr, cnv.gr) {
               standardGeneric("mergeSnpCnv")
           }
)

#' @title mergeSnpCnv
#' @description merge snp and cnv data
#' @rdname mergeSnpCnv
#' @aliases mergeSnpCnv
#' @param .Object the object
#' @param snp.gr SNP GenomicRanges object
#' @param cnv.gr CNV GenomicRanges object
#' @return combined, unique list of genes
setMethod("mergeSnpCnv",
          "RBD",
          function(.Object, snp.gr, cnv.gr) {
              hits = GenomicRanges::findOverlaps(snp.gr, cnv.gr)
              hits.df = as.data.frame(hits)
              snp.df = as.data.frame(snp.gr)
              cnv.df = as.data.frame(cnv.gr)

              colnames(cnv.df) = gsub("start", "cnv.start", colnames(cnv.df))
              colnames(cnv.df) = gsub("end", "cnv.end", colnames(cnv.df))
              cnv.df$seg.id = 1:nrow(cnv.df)

              out = cbind(snp.df[hits.df$queryHits,],
                          cnv.df[hits.df$subjectHits,c("num.mark",
                                                       "seg.mean",
                                                       "seg.id",
                                                       "cnv.start",
                                                       "cnv.end")])
              
              return(out)
          }
)
