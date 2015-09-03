#' @title BTreePlotter
#' @description BTreePlotter
#' @docType package
#' @name BTreePlotter
#' @examples 
#' btreeplotter <- new("BTreePlotter")
utils::globalVariables(c("seg.size", "hds", "lrr", "bubble.size", "R", "HDS", 
                         "Genotype", "symbol.col", "src.seg.size", 
                         "trt.seg.size", "src.lrr", "src.hds", "trt.lrr", 
                         "trt.hds", "src.seqnames", "mutate", "."))

setClass("gg")

BTreePlotter <- setClass(
    "BTreePlotter",

    representation(max.ploidy="numeric",
                   seq.col="character",
                   branch.col="character",
                   branches="gg",
                   max.size="numeric"),

    prototype = list(max.ploidy=6,
                     seq.col=rainbow(22),
                     branch.col="",
                     branches=NULL,
                     max.size=20)
)

setMethod("initialize",
          "BTreePlotter",
          function(.Object, max.ploidy=6, seq.col=NULL,
                   branch.col=NULL, branches=NULL, max.size=20) {

              gg_color_hue <- function(n) {
                  hues = seq(15, 375, length=n+1)
                  hcl(h=hues, l=65, c=100)[1:n]
              }

              if(is.null(seq.col)){
                  seq.col <- rainbow(22)
                  names(seq.col) <- paste0("chr", 1:22)
              }

              if(is.null(branch.col)){
                  branch.col <- gg_color_hue(max.ploidy+1)
                  names(branch.col) <- 0:max.ploidy
              }else{
                  branch.col <- rep(branch.col, length = max.ploidy + 1)
                  names(branch.col) <- 0:max.ploidy
              }

              .Object@max.ploidy = max.ploidy
              .Object@seq.col = seq.col
              .Object@branch.col = branch.col
              .Object@max.size = max.size

              # make the branches
              my.scale <- 5
              xypGrid <- expand.grid(x=0:floor(max.ploidy/2),
                                     y=0:max.ploidy,
                                     p=seq(0,1, by=0.1)) %>% 
                  filter (x <= y & x+y <= max.ploidy ) %>% 
                  unique %>% arrange( -(x+y))

              branches <- xypGrid %>%  
                  dplyr::rowwise() %>% 
                  dplyr::mutate(Genotype=paste0(rep(c("A", "B"), times=c(x,y)), 
                                                collapse=""),
                                R=((x+y)*p + 2*(1-p))/2,
                                HDS=p*(y-x)/4/R) %>% 
                  as.data.frame %>% 
                  mutate(HDS=ifelse(!is.nan(HDS), HDS, 0)) %>% 
                  mutate(symbol.col=ifelse(p==1, "100%", 
                                           ifelse(p==0.5, 
                                                  "50%", 
                                                  "Others" ))) %>% 
                  filter(Genotype != "AB") %>% 
                  mutate(Genotype = ifelse(Genotype == "", "phi", Genotype))

              branches.text <- subset(branches, p ==1 ) %>% 
                  mutate(R= ifelse(Genotype == "phi", R, R+0.1), 
                         HDS = ifelse(Genotype == "phi", HDS+0.03, HDS))

              branches$my.scale <- max.size/2

              # border of the branch + geom_line(aes(cex=p*my.scale+1), alpha=1)
              b0 <- ggplot(branches,
                           aes(x = R, 
                               y=HDS, 
                               group=Genotype, 
                               col=factor(x+y), 
                               fill=symbol.col)) +
                  geom_line(aes(cex=p*my.scale+0.5),
                            alpha=0.1) +
                  geom_point(aes(cex=my.scale*p+4),
                             pch=15, colour = "white") +
                  geom_text(aes(cex=my.scale*p + 0.1,
                                label=round(p*100),
                                angle=ifelse(x!=y, 0, 90)),
                            alpha=1) +
                  geom_point(x=1, y=0, pch=16, size=5, col="gray50") +
                  geom_text(data=branches.text,
                            aes(x=R, y=HDS, label=Genotype),
                            adj=0,
                            parse=TRUE,
                            cex=5,
                            face="bold",
                            family="mono")

              # change theme
              b1 <- b0 + theme(panel.background = element_blank(),
                               axis.title=element_text(size=14,face="bold"),
                               panel.grid.major = element_line(colour="gray90", 
                                                               size=0.3)) +
                  scale_fill_manual(values=c("100%" = "red", 
                                             "50%" = "blue", 
                                             "others" = "white", 
                                             seq.col),
                                    breaks=names(seq.col)) +
                  scale_size_continuous(limits=c(0, max.size), 
                                        range=c(0,10)) +
                  guides(size=FALSE, color=FALSE) + 
                  scale_color_manual(values = c(branch.col, seq.col))

              # change grid line
              bplot <- b1 + scale_y_continuous(breaks = seq(0, 0.5, 0.05)) +
                  scale_x_continuous(breaks = seq(0, max.ploidy/2, 0.5),
                                     limits=c(0, max.ploidy/2+ 0.4))

              .Object@branches <- bplot
              .Object
          }
)

#' @export
#' @docType methods
#' @rdname trackBTree
setGeneric(name="trackBTree",
           def=function(.Object, rbd1, rbd2, is.matched=FALSE,
                        min.srcSize=0.5, min.trtSize=0.1, min.overlap=1e5) {
               standardGeneric("trackBTree")
           }
)

#' @title trackBTree
#' @rdname trackBTree
#' @aliases trackBTree
#' @param .Object the object
#' @param rbd1 rbd one
#' @param rbd2 rbd two
#' @param is.matched is it matched
#' @param min.srcSize min src size
#' @param min.trtSize min trt size
#' @param min.overlap min overlap
#' @return geom_segment location of BTree track
#' @example examples/trackBTree-Ex.R
setMethod(f="trackBTree",
          signature="BTreePlotter",
          definition=function(.Object, rbd1, rbd2, 
                              is.matched=FALSE, min.srcSize=0.5, 
                              min.trtSize=0.1, min.overlap=1e5) {
              if(!is.matched){
                  matched.rbds <- matchSeg(.Object, 
                                           rbd1, 
                                           rbd2, 
                                           min.overlap=min.overlap)
                  rbd1 <- matched.rbds$src
                  rbd2 <- matched.rbds$trt
              }
              dat <- cbind(src=rbd1, trt=rbd2) %>% 
                  filter(src.seg.size>=min.srcSize & trt.seg.size>=min.trtSize)

              if(nrow(dat) == 0) return(NULL)

              out <- geom_segment(data=dat,
                                  aes(x = 2^src.lrr,
                                      y=src.hds,
                                      xend=2^trt.lrr,
                                      yend=trt.hds,
                                      col=src.seqnames,
                                      group=NULL,
                                      fill=src.seqnames),
                                  size=1,
                                  alpha=0.5,
                                  arrow=arrow(length = unit(0.3, "cm")))
              return(out)
          }
)

setGeneric(name="matchSeg",
           def=function(.Object, rbd1, rbd2, min.overlap=1e6) {
               standardGeneric("matchSeg")
           }
)

# find segments
setMethod("matchSeg",
          "BTreePlotter",
          function(.Object, rbd1, rbd2, min.overlap=1e6) {

              rbd1.gr <- rbd1 %>% with(., GRanges(seqnames,
                                                  IRanges(start, end),
                                                  id=seg.id))
              
              rbd2.gr <- rbd2 %>% with(., GRanges(seqnames,
                                                  IRanges(start, end),
                                                  id=seg.id))
              ov <- findOverlaps(rbd1.gr, rbd2.gr, minoverlap = min.overlap) %>% 
                  as.data.frame

              rv <- list(src=rbd1[ov$queryHits, ], trt=rbd2[ov$subjectHits,])
              return(rv)
          }
)

#' @export
#' @docType methods
#' @rdname drawBTree
setGeneric(name="drawBTree",
           def=function(.Object, rbd, size=1) {
               standardGeneric("drawBTree")
           }
)

#' @title drawBTree
#' @rdname drawBTree
#' @aliases drawBTree
#' @param .Object the object
#' @param rbd the rbd object
#' @param size the size
#' @return draw the BTree track
#' @example examples/drawBTree-Ex.R
setMethod("drawBTree",
          "BTreePlotter",
          function(.Object, rbd, size=1) {
              zoom.size <- max(5*size, 
                               .Object@max.size/ max(rbd$seg.size, na.rm=TRUE), 
                               na.rm=TRUE)
              dat <- as.data.frame(rbd) %>% 
                  filter(!is.na(hds)) %>% 
                  mutate(bubble.size = pmin(seg.size * zoom.size, 
                                            .Object@max.size))

              ss <- subset(dat, 2^lrr >= .Object@max.ploidy/2 + 0.4)
              if(nrow(ss)> 0){
                  warning("More ploidy might be suggested: ", 
                          paste(round(2^ss$lrr /2, 1), collapse=", "), "\n")
              }

              btree <- .Object@branches +
                  geom_point(data=dat,
                             aes(x=2^lrr, 
                                 y=hds, 
                                 cex=bubble.size, 
                                 fill=seqnames, 
                                 group=NULL),
                             col="gray50",
                             pch=21,
                             alpha=0.7) +
                  guides(size=FALSE, color=FALSE)

              # add the chr legend
              btree <- btree +
                  guides(fill = guide_legend(override.aes=list(shape = 21,
                                                               alpha=0.7,
                                                               size=4,
                                                               col="black")) ) +
                  theme(legend.key = element_rect(fill=NA), 
                        legend.key.height=grid::unit(0.7,"line")) +
                  ggplot2::labs(fill="Chromosome")
              btree
          }
)

#' @export
#' @docType methods
#' @rdname drawBubbles
setGeneric(name="drawBubbles",
           def=function(.Object, rbd, col=NULL) {
               standardGeneric("drawBubbles")
           }
)

#' @title drawBubbles
#' @rdname drawBubbles
#' @aliases drawBubbles
#' @param .Object the object
#' @param rbd the rbd object
#' @param col the col value
#' @return draw the bubbles on the track
#' @example examples/drawBubbles-Ex.R
setMethod("drawBubbles",
          "BTreePlotter",
          function(.Object, rbd, col="gray80") {
              zoom.size <- max(5, 
                               .Object@max.size/ max(rbd$seg.size, na.rm=TRUE), 
                               na.rm=TRUE)
              
              dat <- as.data.frame(rbd) %>% 
                  filter(!is.na(hds)) %>% 
                  mutate(bubble.size = pmin(seg.size * zoom.size, 
                                            .Object@max.size))

              ss <- subset(dat, 2^lrr >= .Object@max.ploidy/2 + 0.4)
              if(nrow(ss)> 0){
                  warning("More ploidy might be suggested: ", 
                          paste(round(2^ss$lrr /2, 1), collapse=", "), "\n")
              }

              if(is.null(col)){
                  bubbles <- geom_point(data=dat,
                                        aes(x=2^lrr,
                                            y=hds,
                                            cex=bubble.size,
                                            fill=seqnames,
                                            group=NULL),
                                        col=col,
                                        pch=21,
                                        alpha=0.7) +
                      guides(size=FALSE, color=FALSE, fill=FALSE)
              }else{
                  bubbles <- geom_point(data=dat,
                                        aes(x=2^lrr, 
                                            y=hds, 
                                            cex=bubble.size, 
                                            group=NULL),
                                        col="gray50",
                                        fill=col,
                                        pch=21,
                                        alpha=0.7)
              }
              return(bubbles)
          }
)

#' @export
#' @docType methods
#' @rdname drawFeatures
setGeneric(name="drawFeatures",
           def=function(.Object, rbd, col=NULL) {
               standardGeneric("drawFeatures")
           }
)

#' @title drawFeatures
#' @rdname drawFeatures
#' @aliases drawFeatures
#' @param .Object the object
#' @param data additional annotation to plot
#' @param col the col value
#' @return draw the annotation on the track
#' @example examples/drawFeatures-Ex.R
setMethod("drawFeatures",
          "BTreePlotter",
          function(.Object, rbd, col="black") {
              
              dat <- as.data.frame(rbd) %>% filter(!is.na(hds)) 
              
              features <- geom_point(data=dat,
                                     aes(x=(2^lrr) + .15,
                                         y=hds, 
                                         cex=15, 
                                         group=NULL),
                                     col="black",
                                     fill=col,
                                     pch="-")
              
              return(features)
          }
)
