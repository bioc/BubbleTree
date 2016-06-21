#' @title BTreePredictor
#' @description BTreePredictor
#' @docType package
#' @name BTreePredictor
#' @examples 
#' btreepredictor <- new("BTreePredictor")
utils::globalVariables(c("seg.size", "seg.id", "x", "y", "hds", "lrr", "R", 
                         "HDS", "Genotype", "symbol.col", "cls", ".p.dif", 
                         "deviation", "prev", "."))

BTreePredictor <- setClass(
    "BTreePredictor",

    representation(config="list",
                   rbd="data.frame",
                   rbd.adj="data.frame",
                   result="list"),

    prototype(config=NULL,
              rbd=NULL,
              rbd.adj=NULL,
              result=NULL)
)

setMethod("initialize",
          "BTreePredictor",
          function(.Object, config=NULL, rbd=NULL, rbd.adj=NULL, result=NULL,
                   max.ploidy=6, prev.grid=seq(0.2,1, by=0.01)) {

              if(is.null(config)) {
                  config <- list()
                  config$xypGrid = expand.grid(x=0:floor(max.ploidy/2),
                                               y=0:max.ploidy,
                                               p=prev.grid) %>% 
                      filter (x <= y & x+y <= max.ploidy) %>% 
                      filter( ! (x==1 & y==1) ) %>% 
                      unique

                  config <- c(config,
                              list(max.distL = 0.06,
                                   max.distR = 0.03,
                                   best.tol = 0.005,
                                   root.hdsInt=0.05,
                                   root.rInt=0.3,
                                   min.segSize = 0.5,
                                   min.prev=0.15,
                                   cutree.h=0.2,
                                   high.ploidy=TRUE,
                                   high.purity=TRUE,
                                   total.mark=NA,
                                   cnv.gr=NULL,
                                   max.hds.sd=0.3,
                                   verbose=TRUE,
                                   min.het.cnt = 20,
                                   lowest.purity=0.4)
                  )
              }

              .Object@config <- config
              .Object@rbd <- data.frame()
              .Object@rbd.adj <- data.frame()
              .Object@result <- list()
              
              if( !is.null(rbd) ) {
                  .Object@rbd <- rbd
              }
              
              .Object
          }
)

#' @export
#' @docType methods
#' @rdname loadRBD
setGeneric(name="loadRBD",
           def=function(.Object, rbd, total.mark=NA) {
               standardGeneric("loadRBD")
           }
)

#' @title loadRBD
#' @description load the RBD data
#' @rdname loadRBD
#' @aliases loadRBD
#' @param .Object the object
#' @param rbd rbd object
#' @param total.mark total mark
#' @return .Object populated with the RBD list with updated segment size
#' @example examples/loadRBD-Ex.R
setMethod("loadRBD",
          "BTreePredictor",
          function(.Object, rbd, total.mark=NA) {
              rbd <- as.data.frame(rbd@elementMetadata)
              cfg <- .Object@config
              if(is.na(total.mark)){
                  total.mark <- sum(rbd$num.mark, na.rm=TRUE)
              }

              rbd$seg.size <- rbd$num.mark / total.mark * 100
              .Object@rbd <- rbd

              .Object@rbd.adj <- data.frame()
              .Object@result <- list()
              .Object
          }
)

#' @export
#' @docType methods
#' @rdname btpredict
setGeneric(name="btpredict",
           signature=c(".Object" ),
           def=function(.Object) {
               standardGeneric("btpredict")
           }
)

#' @title btpredict
#' @description btpredict
#' @rdname btpredict
#' @aliases btpredict
#' @param .Object the object
#' @return .Object populated with the predictions
#' @example examples/btpredict-Ex.R
setMethod("btpredict",
          "BTreePredictor",
          function(.Object) {
              .Object <- getPloidyAdj(.Object)
              
              # find the prevalences of the (sub)clones
              .Object <- findClones(.Object, useABB=FALSE)
              
              if( is.na(.Object@config$high.ploidy)) {
                  .Object@config$high.ploidy = TRUE
              }
              
              ###############################
              # Extra calculation for AABB
              ###############################
              round=1
              if(.Object@result$ploidy.adj["ploidy"] == 4 &
                 .Object@config$high.ploidy) {
                  # the tumor ploidy could be AABB or AAABBB, etc.
                  # if purity is less than 1/3, there is a chance of AAABBB
                  # or higher but we only consider AABB here
                  p <- .Object@result$prev
                  p <- p[p<0.5]
                  if(length(p) == 0 || p == 0 || is.na(p)){
                      # highPloidy estimation could be wrong
                      .Object@result$ploidy.adj["ploidy"] <- 2
                  }else{
                      new.adj <- 1/(1-max(p))
                      .Object@result$ploidy.adj["adj"] <- 
                          .Object@result$ploidy.adj["adj"] * new.adj
                      .Object@rbd.adj$lrr <- .Object@rbd$lrr +
                          log2(.Object@result$ploidy.adj["adj"])
                      round <- 2
                      .Object <- findClones(.Object, useABB=TRUE)
                  }
              }
              
              ###############################
              # last, find the best match of all segments including
              # indistinguiable branches
              ###############################
              pp <- sort(.Object@result$prev, decreasing=TRUE)
              
              if(.Object@result$ploidy.adj["ploidy"] == 3 &
                 pp[1] < .Object@result$ploidy.adj["purity"] -
                 .Object@config$cutree.h) {
                  pp <- c(.Object@result$ploidy.adj["purity"], pp)
              }
              
              if(!is.na(pp[1]) & pp[1] > .Object@config$min.prev) {
                  # most sublones below 30% could be artifacts as
                  # hds is not be exactly
                  pp <- pp[pp> .Object@config$min.prev]
              }
              pred0 <- c(pp)
              
              predFit <- new("BTreePredictor", 
                             rbd=.Object@rbd, 
                             rbd.adj=.Object@rbd.adj, 
                             prev.grid=pred0, 
                             max.ploidy=10)

              predFit@config$best.tol <- 0.02
              predFit@config$max.distL <- 999
              predFit@config$max.distR <- 999
              predFit@config$xypGrid <- rbind(predFit@config$xypGrid, c(1,1,1))
              
              .Object@result$dist <- adply(.Object@rbd.adj, 
                                           1, 
                                           function(df) 
                                               findBestXYP(predFit, 
                                                           df[1, "lrr"], 
                                                           df[1, "hds"])) %>% 
                  plyr::ddply(.(seg.id), 
                              function(df) df %>% 
                                  filter(abs(x-1) + 
                                             abs(y-1) == min(abs(x-1) + 
                                                                 abs(y-1))) %>% 
                                  filter(dist == min(dist, na.rm=TRUE)) )

              .Object@result$prev <- pp

              deviation <- with(subset(.Object@result$dist, seg.size > 0.5),
                                sum(seg.size * dist/100, na.rm=TRUE))
              .Object@result$deviation <- deviation

              return(.Object)
          }
)

setGeneric(name="getPloidyAdj",
           def=function(.Object) {
               standardGeneric("getPloidyAdj")
           }
)

setMethod("getPloidyAdj",
          "BTreePredictor",
          function(.Object) {
              cfg <- .Object@config
              rbd <- .Object@rbd

              if( is.na(cfg$high.ploidy) ) {
                  cfg$high.ploidy = TRUE
              }

              total.size <- sum(rbd$seg.size, na.rm=TRUE)
              centralSegInterval = 0.1
              uc.cov.cutoff=0.5
              lowHds.cutoff = 0.15
              uc.hds.cutoff=0.07
              highPloidy=cfg$high.ploidy

              all.seg.mean <- with(.Object@rbd,
                                   limma::weighted.median(lrr,
                                                          seg.size,
                                                          na.rm=TRUE))

              if(abs(all.seg.mean) > 0.3){
                  # special abnormal case like 1500
                  central.segs <- subset(rbd, abs(lrr) < centralSegInterval)
                  all.seg.mean <- with(central.segs,
                                       limma::weighted.median(lrr,
                                                              seg.size,
                                                              na.rm=TRUE))
              }else{
                  central.segs <- subset(rbd,
                                         abs(lrr - all.seg.mean) < 
                                                        centralSegInterval)
              }
              
              uc.segs <- subset(central.segs, hds > uc.hds.cutoff)

              uc.cov <- sum(uc.segs$seg.size, na.rm=TRUE) /
                            sum(central.segs$seg.size, na.rm=TRUE)
            
              low.hds <- with(rbd,
                              sum(seg.size[hds < lowHds.cutoff], na.rm=TRUE)) /
                              total.size
              
              info <- c(upperCentral.cov = uc.cov, lowHDS.cov = low.hds)

              ploidy <- 2
              adj <- 2^all.seg.mean
              uc.hds <- 0
              if(nrow(uc.segs) >0) {
                  uc.hds <- limma::weighted.median(uc.segs$hds,
                                                   uc.segs$seg.size,
                                                   na.rm=TRUE)
              }

              if( is.na(cfg$high.ploidy)) {
                  cfg$high.ploidy=TRUE
              }

              if(uc.cov > 0.5 & uc.hds < 1/6 & cfg$high.ploidy){
                  ploidy <- 3
                  uc.lrr <- with(uc.segs, limma::weighted.median(lrr, seg.size))
                  purity <- 4*uc.hds/(1-2*uc.hds)
                  adj <- (1+purity/2) /2^uc.lrr

                  # also check the new central stages to make the 
                  # fine adjustment
                  central.segs.new <- subset(rbd,
                                             abs(lrr + log2(adj)-all.seg.mean) < 
                                                 centralSegInterval)
                  if(nrow(central.segs.new) > 0 ){
                      adj <-  with(central.segs.new,
                                   2^(-limma::weighted.median(lrr, seg.size)))
                  }
              }else{
                  # do not estimate the purity here
                  purity <- NA
                  lrr <- with(central.segs, weighted.mean(lrr, seg.size))
                  adj <- 1 / 2^lrr

                  ploidy <- 2
                  if(low.hds >0.95 & highPloidy ){
                      #and no seg in the upper right
                      upright.segs <- subset(rbd, hds > 0.1 & 2^lrr > 1.25)
                      if( sum(upright.segs$seg.size, na.rm=TRUE)  < 1){
                          ploidy <- 4
                      }
                  }
              }

              .Object@result$ploidy.adj <- c(info,
                                             adj=adj,
                                             ploidy=as.integer(ploidy),
                                             purity=purity)

              rbd.adj <- .Object@rbd
              rbd.adj$lrr <- rbd.adj$lrr + log2(adj)
              .Object@rbd.adj <- rbd.adj
              return(.Object)
          }
)

setGeneric(name="findClones",
           def=function(.Object, useABB=FALSE) {
               standardGeneric("findClones")
           }
)

setMethod("findClones",
          "BTreePredictor",
          function(.Object, useABB=FALSE) {
              cfg <- .Object@config
              rbd.adj <- .Object@rbd.adj

              # filter some ploidy states not informative to identify
              # clones and also save the times
              .Object@config$xypGrid <- .Object@config$xypGrid %>% 
                  filter ( !((x==0 & y >= 5) | (x>1 & y-x == 1)))
              
              rbd.flt <- subset(rbd.adj,
                                seg.size > cfg$min.segSize & 
                                    !(hds < cfg$root.hdsInt & abs(lrr) < 
                                          cfg$root.rInt))

              # with the distance
              rbd.dist <- plyr::adply(rbd.flt,
                                      1,
                                      function(df) findBestXYP(.Object, 
                                                               df[1, "lrr"], 
                                                               df[1, "hds"]))

              if(nrow(rbd.dist) == 0) {
                  # no prediction
                  .Object@result$prev <- NA
                  return(.Object)
              }

              # now predict the purity and
              #
              rbd.dist.flt <- subset(rbd.dist, x!=1 & x!=y & p > cfg$min.prev )
              is.ambi <- FALSE

              if( is.na(cfg$high.purity)) {
                  cfg$high.purity = TRUE
              }

              if(nrow(rbd.dist.flt) == 0 || useABB || 
                 (cfg$high.purity && 
                  max(rbd.dist.flt$p, na.rm=TRUE) < cfg$lowest.purity)) {
                  # check the unidentifiable states
                  # and only use the one with lowest y
                  rbd.dist.flt <- subset(rbd.dist,
                                         p > cfg$min.prev) %>% 
                      plyr::ddply(.(seg.id), 
                                  function(df) df %>% 
                                      filter(dist<min(dist)+cfg$best.tol) %>% 
                                      filter(y == min(y) ) )
                  is.ambi <- TRUE
              }

              if(nrow(rbd.dist.flt) == 0){
                  .Object@result$prev <- NA #data.frame(cls=1, prev=NA)
                  return(.Object)
              }else{
                  if(nrow(rbd.dist.flt)==1){
                      prev <- data.frame(cls=1, prev=rbd.dist.flt$p)
                  }else{
                      rbd.dist.flt$cls <- cutree(hclust(dist(rbd.dist.flt$p)),
                                                 h=cfg$cutree.h)

                      # drop those cluster with few segment
                      if(is.ambi){
                          abb.segs <- rbd.dist.flt %>% group_by(cls) %>% 
                              filter(all(x ==1))
                          
                          if(nrow(abb.segs) > 0) {
                              # if there are some ABB-only segments, try to 
                              # drop them if the coverage is low
                              abb.segs  <- abb.segs %>% 
                                  dplyr::summarise(cov=sum(seg.size, 
                                                           na.rm=TRUE)) %>% 
                                  filter(cov < 0.5)

                              if(length(abb.segs$cls) > 0) {
                                  rbd.dist.flt <- subset(rbd.dist.flt,
                                                         !cls %in% abb.segs$cls)
                              }
                          }
                      }
                      prev <- group_by(rbd.dist.flt, cls) %>% 
                          dplyr::summarise(prev=weighted.mean(p, seg.size)) %>% 
                          as.data.frame
                  }
              }

              # new
              if(!is.ambi){
                  # make sure all the major ABB states covered by at least 
                  # one of prev add the one with the smallest likely prev 
                  # otherwise at least 1%
                  rbd.dist.ambi <- subset(rbd.dist, 
                                          (x==1 | x==y) & p > 
                                              max(1,cfg$min.prev) )

                  if(nrow(rbd.dist.ambi) > 0) {
                      # prev on the ABB  is more sensitive
                      rbd.dist.ambi.out <- rbd.dist.ambi %>% 
                          mutate(.p.dif = 
                                     sapply(p, 
                                            function(x) min (abs(x-prev$prev), 
                                                             na.rm=TRUE))) %>% 
                          group_by(seg.id) %>% 
                          filter( all( .p.dif > 2*cfg$cutree.h))
                      
                      if(nrow(rbd.dist.ambi.out)>0){

                          rbd.dist.ambi.out <- rbd.dist.ambi.out %>% 
                              group_by(seg.id) %>% 
                              filter(x == min(x))
                          
                          if(nrow(rbd.dist.ambi.out) == 1){
                              prev <- rbind(prev, 
                                            c(nrow(prev)+1, 
                                              rbd.dist.ambi.out$p))
                          }else {
                              rbd.dist.ambi.out$cls <- 
                                  cutree(hclust(dist(rbd.dist.ambi.out$p)),
                                         h=cfg$cutree.h)
                              prev <- group_by(rbd.dist.ambi.out, cls) %>% 
                                  dplyr::summarise(
                                      prev = weighted.mean(p, seg.size)) %>% 
                                  as.data.frame %>% 
                                  rbind(prev, .)
                          }
                      }
                  } # end of
              } # end of is.ambi

              .Object@result$prev <- prev$prev
              .Object@result$dist <- rbd.dist

              .Object
          }
)

setGeneric(name="findBestXYP",
           signature=c(".Object" ),
           function(.Object, lrr, hds) {
               standardGeneric("findBestXYP")
           }
)

setMethod("findBestXYP",
          "BTreePredictor",
          function(.Object, lrr, hds) {
              cfg <- .Object@config

              WW = 5 # weight to reduce the effect R sccore
              out <- plyr::ddply(cfg$xypGrid, .(x,y,p), function(df){
                  x <- df[1,"x"]
                  y <- df[1,"y"]
                  p <- df[1,"p"]
                  tau <- (x+y)*p + 2*(1-p)
                  lrr.p <- log2(tau) - 1
                  hds.p <- ifelse(y==x, 0, p*(y-x)/2/tau)
                  dist <- sqrt( (2^lrr - tau/2)^2/WW + (hds - hds.p)^2)
                  return(c(dist=dist, lrr.pred = lrr.p, hds.pred = hds.p))
              })

              #only keep the best one for each x y combo
              out <- subset(out,
                            (lrr < 0.1 & dist < cfg$max.distL) | 
                                (lrr>=0.1 & dist < cfg$max.distR))

              na.out <- data.frame(x=NA, y=NA, p=NA, dist=NA)

              if(nrow(out) == 0) return(na.out)

              out <- out %>% group_by(x, y) %>% filter (dist == min(dist))

              if(nrow(out) >= 1){
                  out <- subset(out, dist <= min(dist) + cfg$best.tol)
              }else{
                  out <- na.out
              }
              return(out)
          }
)

#' @export
#' @docType methods
#' @rdname info
setGeneric(name="info",
           def=function(.Object) {
               standardGeneric("info")
           }
)

#' @title info
#' @description info
#' @rdname info
#' @aliases info
#' @param .Object the object
#' @return print out info of prediction data
#' @example examples/info-Ex.R
setMethod("info",
          "BTreePredictor",
          function(.Object) {
              if(is.null(.Object@result$dist)){
                  info <- "To be predicted!"
                  return(info)
              }

              purity <- .Object@result$prev[1]
              adj <- .Object@result$ploidy.adj["adj"]
              
              # when purity is low the calculation result is not reliable
              ploidy <- (2*adj -2)/purity + 2
              info <- with(.Object@result,
                           sprintf(
                               "Purity: %s; Ploidy: %3.1f; Deviation: %4.2f",
                               paste(round(prev,2), collapse=", "),
                               round(ploidy,2), deviation))
              return(info)
          }
)

