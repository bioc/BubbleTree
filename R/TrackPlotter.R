#' @title TrackPlotter
#' @description TrackPlotter
#' @docType package
#' @name TrackPlotter
#' @examples
#' trackplotter <- new("TrackPlotter")
utils::globalVariables(c("p", "x", "y", "chrom", "chromStart", "chromEnd", 
                         ".start", ".end", "r.pred", "hds.pred", "x.adj", 
                         "y.adj", "brewer.pal", "metadata", "seqlevels<-", 
                         "mutate"))




gr2 <- setClass("GenomicRanges")

TrackPlotter <- setClass(
    "TrackPlotter",

    representation(result.dat="data.frame"),

    prototype = list(result.dat=NULL)
)

setMethod("initialize",
          "TrackPlotter",
          function(.Object, result.dat=NULL) {
              .Object@result.dat <- data.frame()
              .Object
          }
)

#' @export
#' @docType methods
#' @rdname xyTrack
setGeneric(name="xyTrack",
           def=function(.Object, result.dat, gr2, min.prev=0.15, ymax=4.3) {
               standardGeneric("xyTrack")
           }
)

#' @title xyTrack
#' @rdname xyTrack
#' @aliases xyTrack
#' @param .Object the object
#' @param result.dat result dataframe
#' @param gr2 gr2 object
#' @param min.prev previous min
#' @param ymax the max y
#' @return the highlighted xy track
#' @example examples/xyTrack-Ex.R
setMethod("xyTrack",
          "TrackPlotter",
          function(.Object, result.dat, gr2, min.prev=0.15, ymax=4.3){
              dat <-result.dat$dist
              main.prev <- result.dat$prev[1]
              subclone.prev <- result.dat$prev[2]

              cols <- brewer.pal(7,"Set1")
              baseline.adj <- 0.1

              # convert 1,1,1 to 1,1, prev
              dat.2 = dat %>%
                  mutate(x=ifelse(p < min.prev, 1, x), 
                         y=ifelse(p < min.prev, 1, y)) %>%
                  mutate(y.adj = y+baseline.adj, x.adj = x-baseline.adj) %>%
                  mutate (p = ifelse(x ==1 & y==1 & p ==1, main.prev, p))
    
              gr1 <- dat2gr(dat.2) %>% keep24chr(no.sex.chr=TRUE)
              gr1.t <- gr1 %>% transformToGenome(space.skip=0.01)

              gr2 <- with(gr2, 
                          GRanges(chrom, IRanges(chromStart, chromEnd))) %>% 
                  keep24chr(no.sex.chr=TRUE)
              
              gr2.t <- gr2 %>% transformToGenome(space.skip=0.01)

              # only label chr1 to 22
              ss2 <- mold(gr2.t) %>% filter(seqnames %in% paste0('chr', 1:22)) 

              my.dat <- mold(gr1.t)
              p1 <- ggplot(my.dat, 
                            aes(x=.start, xend=.end, y=y.adj, yend=y.adj)) +
                            ylab("Copy number") +
                            xlab("Genomic location") +
                            geom_hline(yintercept = 1, 
                                       colour = "gray", 
                                       lty=1, 
                                       alpha=0.2, 
                                       size = 0.3) +
                            geom_segment(col=cols[1], size=1, na.rm=TRUE)  +
                            geom_segment(data=my.dat, 
                                         aes(x=.start, 
                                             xend=.end, 
                                             y=x.adj, 
                                             yend=x.adj), 
                                         col=cols[2], 
                                         size=1, 
                                         inherit.aes=FALSE,
                                         na.rm=TRUE) +
                            #scale_y_continuous(breaks = seq(0,9,1)) +
                            theme(panel.background = element_blank()) +
                            geom_hline(yintercept = 1, 
                                       colour = "gray80", 
                                       lty=1, 
                                       size = 0.5)

              p2 <- p1 + scale_x_continuous(breaks = (ss2$.start + ss2$.end)/2,
                                            labels=as.character(ss2$seqnames)) +
                  theme(panel.background = element_blank(),
                        panel.grid.major = element_line(colour="gray90", 
                                                        size=0.3),
                        panel.grid.minor = element_blank(),
                        axis.text.x = element_text(angle = 45, 
                                                   hjust = 1, 
                                                   size = 10)) +
                  ylim(c(-0.3,ymax))

              if(!is.na(subclone.prev)){
                  sp <- result.dat$prev[-1]
                  dd <- subset(my.dat, p <= subclone.prev)
                  cc <- brewer.pal(8,"Set2")[1:length(sp)]
                  fill <- plyr::mapvalues(dd$p, from=sp, to=cc)
                  p3 <- p2 + geom_rect(data=dd, 
                                       aes(xmin=.start, 
                                           xmax=.end, 
                                           ymin=-Inf, 
                                           ymax=Inf), 
                                       fill=fill, 
                                       alpha=0.2, 
                                       inherit.aes = FALSE)

              }else{
                  p3 <- p2
              }

              # high light the subclone
              return(p3)
          }
)

#' @export
#' @docType methods
#' @rdname bafTrack
setGeneric(name="bafTrack",
           def=function(.Object, result.dat, gr2, 
                        somatic.gr=NULL, min.prev=0.15, cex=1.2) {
               standardGeneric("bafTrack")
           }
)

#' @title bafTrack
#' @rdname bafTrack
#' @aliases bafTrack
#' @param .Object the object
#' @param result.dat the result dataframe
#' @param gr2 the gr2 object
#' @param somatic.gr somatic gr object annotation
#' @param min.prev previous min
#' @param cex the cex
#' @return the highlighted BAF track
#' @example examples/bafTrack-Ex.R
setMethod("bafTrack",
          "TrackPlotter",
          function(.Object, result.dat, gr2, 
                   somatic.gr=NULL, min.prev=0.15, cex=1.2){
              dat <-  result.dat$dist
              main.prev <- result.dat$prev[1]
              subclone.prev <- result.dat$prev[2]

              # assume gr1 and gr2 has similar seqinfo
              #seqlevels(gr2) <-  seqlevels(gr1)
              #seqinfo(gr2) <- seqinfo(gr1)
              baseline.adj <- 0.1

              cols <- brewer.pal(8,"Dark2")
              # convert 1,1,1 to 1,1, prev
              dat.2 = dat %>% 
                  mutate(x = ifelse(p < min.prev, 1, x), 
                         y=ifelse(p < min.prev, 1, y)) %>% 
                  mutate(y.adj = y+baseline.adj, 
                         x.adj = x-baseline.adj) %>% 
                  mutate (p = ifelse(x==1 & y==1 & p > main.prev, main.prev, p))

              gr1 <- dat2gr(dat.2) %>% keep24chr(no.sex.chr=TRUE)
              gr1.t <- gr1 %>% transformToGenome(space.skip=0.01)

              gr2 <- with(gr2, GRanges(chrom, 
                                       IRanges(chromStart, chromEnd))) %>% 
                  keep24chr(no.sex.chr=TRUE)
              
              gr2.t <- gr2 %>% transformToGenome(space.skip=0.01)

              # only label chr1 to 22
              ss2 <- mold(gr2.t) %>% filter(seqnames %in% paste0('chr', 1:22)) 

              my.dat <- mold(gr1.t)
              my.dat$main.prev <- main.prev

              p0 <- ggplot(inherit.aes=FALSE)
              if(!is.null(somatic.gr)){
                  somatic.dat <- somatic.gr %>% keep24chr(no.sex.chr=TRUE)  %>% 
                      transformToGenome(space.skip=0.01) %>% mold
                  p0 <- p0 + geom_point(data=somatic.dat, 
                                        aes(x=.start, y=score), 
                                        inherit.aes=FALSE, 
                                        pch=21, 
                                        cex=cex, 
                                        col="black", 
                                        fill=cols[7], 
                                        alpha=0.5)
              }

              p1 <-  p0 + geom_segment(data=my.dat, 
                                       aes(x=.start, 
                                           xend=.end, 
                                           y = (p*x)/ ( (x+y)*p + 2*(1-p)), 
                                           yend=(p*x)/( (x+y)*p + 2*(1-p))), 
                                       col=cols[4], size=1, lty=1, na.rm=TRUE) + 
                  geom_segment(data=my.dat, 
                               aes(x=.start, 
                                   xend=.end, 
                                   y = (main.prev-p+p*x)/ ( (x+y)*p + 2*(1-p)), 
                                   yend=(main.prev-p+p*x)/ ((x+y)*p + 2*(1-p))), 
                               col=cols[3], size=1, lty=1, na.rm=TRUE) +  
                  geom_segment(data=my.dat, 
                               aes(x=.start, 
                                   xend=.end, 
                                   y = y*p/ ( (x+y)*p + 2*(1-p)), 
                                   yend=y*p/ ( (x+y)*p + 2*(1-p))), 
                               col=cols[2], 
                               size=1, 
                               lty=1, na.rm=TRUE) + 
                  geom_segment(data=my.dat, 
                               aes(x=.start, 
                                   xend=.end, 
                                   y = (y*p + main.prev-p)/((x+y)*p + 2*(1-p)), 
                                   yend=(y*p + main.prev-p)/((x+y)*p+2*(1-p))), 
                               col=cols[1], 
                               size=1, 
                               lty=1, na.rm=TRUE) + 
                  ylim(c(0,1)) + 
                  ylab("Max(BAF)")

              p2 <- p1 + scale_x_continuous(breaks = (ss2$.start + ss2$.end)/2,
                                            labels=as.character(ss2$seqnames)) + 
                  theme(panel.background = element_blank(), 
                        panel.grid.major = element_line(colour="gray90", 
                                                        size=0.3), 
                        panel.grid.minor = element_blank(), 
                        axis.text.x = element_text(angle=45, hjust=1, size=10))

              if(!is.na(subclone.prev)){
                  sp <- result.dat$prev[-1]
                  dd <- subset(my.dat, p <= subclone.prev)
                  cc <- brewer.pal(8,"Set2")[1:length(sp)]
                  fill <- plyr::mapvalues(dd$p, from=sp, to=cc)
                  p3 <- p2 + geom_rect(data=dd, 
                                       aes(xmin=.start, 
                                           xmax=.end, 
                                           ymin=-Inf, 
                                           ymax=Inf), 
                                       fill=fill, 
                                       alpha=0.2, 
                                       inherit.aes = FALSE)

              }else{
                  p3 <- p2
              }

              return(p3)
          }
)

#' @export
#' @docType methods
#' @rdname heteroLociTrack
setGeneric(name="heteroLociTrack",
           def=function(.Object, result.dat, gr2, 
                        hetero.gr=NULL, min.prev=0.15, ymax=4.3, cex=0.5) {
               standardGeneric("heteroLociTrack")
           }
)

#' @title heteroLociTrack
#' @rdname heteroLociTrack
#' @aliases heteroLociTrack
#' @param .Object the object
#' @param result.dat the results
#' @param gr2 the gr2 object
#' @param hetero.gr hetero annotation
#' @param min.prev previous min
#' @param ymax max y
#' @param cex the cex
#' @return the highlightted heterozygosity track
#' @example examples/heteroLociTrack-Ex.R
setMethod("heteroLociTrack",
          "TrackPlotter",
          function(.Object, result.dat, gr2, 
                   hetero.gr=NULL, min.prev=0.15, ymax=4.3, cex=0.5){
              dat <-  result.dat$dist
              main.prev <- result.dat$prev[1]
              subclone.prev <- result.dat$prev[2]

              # assume gr1 and gr2 has similar seqinfo
              #seqlevels(gr2) <-  seqlevels(gr1)
              #seqinfo(gr2) <- seqinfo(gr1)
              baseline.adj <- 0.1

              cols <- brewer.pal(8,"Dark2")
              col.p <- brewer.pal(8,"Paired")

              # convert 1,1,1 to 1,1, prev
              dat.2 <- dat %>% mutate(x = ifelse(p < min.prev, 1, x), 
                                      y=ifelse(p < min.prev, 1, y)) %>% 
                  mutate(y.adj = y+baseline.adj, x.adj = x-baseline.adj) %>% 
                  mutate (p=ifelse(x ==1 & y==1 & p > main.prev, main.prev, p))

              # calculate the predicted hds (this step can be omitted if 
              # findBestXYP has been updated)
              dat.2 <- dat.2 %>% 
                  mutate(hds.pred = p*(y-x)/2/((x+y)*p + 2*(1-p)) )

              gr1 <- dat2gr(dat.2) %>% keep24chr(no.sex.chr=TRUE)
              gr1.t <- gr1 %>% transformToGenome(space.skip=0.01)

              gr2 <- with(gr2, GRanges(chrom, 
                                       IRanges(chromStart, chromEnd))) %>% 
                  keep24chr(no.sex.chr=TRUE)
              gr2.t <- gr2 %>% transformToGenome(space.skip=0.01)

              # only label chr1 to 22
              ss2 <- mold(gr2.t) %>% filter(seqnames %in% paste0('chr', 1:22)) 

              my.dat <- mold(gr1.t)
              my.dat$main.prev <- main.prev
              # if(is.na(main.prev)) return(ggplot() + geom_blank())

              p0 <- ggplot(inherit.aes=FALSE)
              if(!is.null(hetero.gr)){
                  hetero.dat <- hetero.gr %>% keep24chr(no.sex.chr=TRUE) %>% 
                      transformToGenome(space.skip=0.01) %>% mold
                  p0 <- p0 + geom_point(data=hetero.dat, 
                                        aes(x=.start, y=score), 
                                        inherit.aes=FALSE, 
                                        pch=21, 
                                        cex=cex, 
                                        col=col.p[2], 
                                        fill=col.p[1], 
                                        alpha=0.2)
              }

              p1 <- p0 +  geom_segment(data=my.dat,
                                       aes(x=.start, 
                                           xend=.end, 
                                           y=0.5 - hds.pred, 
                                           yend=0.5 - hds.pred), 
                                       col=cols[2], 
                                       size=2, 
                                       lty=1, na.rm=TRUE) + 
                  geom_segment(data=my.dat, 
                               aes(x=.start, 
                                   xend=.end, 
                                   y = 0.5 + hds.pred, 
                                   yend=0.5 + hds.pred), 
                               col=cols[2], 
                               size=2, 
                               lty=1, na.rm=TRUE)  + 
                  ylim(c(0,1)) + 
                  ylab("BAF")

              p2 <- p1 + scale_x_continuous(breaks = (ss2$.start + ss2$.end)/2,
                                            labels=as.character(ss2$seqnames)) + 
                  theme(panel.background = element_blank(), 
                        panel.grid.major = element_line(colour="gray90", 
                                                        size=0.3), 
                        panel.grid.minor = element_blank(), 
                        axis.text.x = element_text(angle=45, hjust=1, size=10))

              if(!is.na(subclone.prev)){
                  sp <- result.dat$prev[-1]
                  dd <- subset(my.dat, p <= subclone.prev)
                  cc <- brewer.pal(8,"Set2")[1:length(sp)]
                  fill <- plyr::mapvalues(dd$p, from=sp, to=cc)
                  p3 <- p2 + geom_rect(data=dd, 
                                       aes(xmin=.start, 
                                           xmax=.end, 
                                           ymin=-Inf, 
                                           ymax=Inf), 
                                       fill=fill, 
                                       alpha=0.2, 
                                       inherit.aes = FALSE)

              }else{
                  p3 <- p2
              }

              return(p3)
          }

)

#' @export
#' @docType methods
#' @rdname RscoreTrack
setGeneric(name="RscoreTrack",
           def=function(.Object, result.dat, gr2, 
                        cnv.gr=NULL, min.prev=0.15, ymax=3, cex=1.5) {
               standardGeneric("RscoreTrack")
           }
)

#' @title RscoreTrack
#' @rdname RscoreTrack
#' @aliases RscoreTrack
#' @param .Object the object
#' @param result.dat the results
#' @param gr2 the gr2 object
#' @param cnv.gr cnv annotation
#' @param min.prev previous min
#' @param ymax max y
#' @param cex the cex
#' @return the highlighted RScore track
#' @example examples/RscoreTrack-Ex.R
setMethod("RscoreTrack",
          "TrackPlotter",
          function(.Object, result.dat, gr2, 
                   cnv.gr=NULL, min.prev=0.15, ymax=3, cex=1.5){
              dat <-  result.dat$dist
              main.prev <- result.dat$prev[1]
              subclone.prev <- result.dat$prev[2]

              # assume gr1 and gr2 has similar seqinfo
              #seqlevels(gr2) <-  seqlevels(gr1)
              #seqinfo(gr2) <- seqinfo(gr1)
              baseline.adj <- 0.1

              cols <- brewer.pal(8,"Dark2")
              # convert 1,1,1 to 1,1, prev
              dat.2 <- dat %>% mutate(x = ifelse(p < min.prev, 1, x), 
                                      y=ifelse(p < min.prev, 1, y)) %>% 
                  mutate(y.adj = y+baseline.adj, x.adj = x-baseline.adj) %>% 
                  mutate (p = ifelse(x ==1 & y==1 & p > main.prev,main.prev, p))

              # calculate the predicted hds (this step can be omitted if 
              # findBestXYP has been updated)
              dat.2 <- dat.2 %>% 
                  mutate(r.pred=((x+y)*p/2 +(1-p))/result.dat$ploidy.adj["adj"])

              gr1 <- dat2gr(dat.2) %>% keep24chr(no.sex.chr=TRUE)

              gr1.t <- gr1 %>% transformToGenome(space.skip=0.01)

              gr2 <- with(gr2, GRanges(chrom, 
                                       IRanges(chromStart, chromEnd))) %>% 
                  keep24chr(no.sex.chr=TRUE)
              
              gr2.t <- gr2 %>% transformToGenome(space.skip=0.01)

              # only label chr1 to 22
              ss2 <- mold(gr2.t) %>% filter(seqnames %in% paste0('chr', 1:22))

              my.dat <- mold(gr1.t)
              my.dat$main.prev <- main.prev

              p0 <- ggplot(inherit.aes=FALSE)
              if(!is.null(cnv.gr)){
                  cnv.gr <- cnv.gr %>% keep24chr(no.sex.chr=TRUE) %>% 
                            transformToGenome(space.skip=0.01) %>% mold
                  p0 <- p0 + geom_segment(data=cnv.gr, 
                                          aes(x=.start, 
                                              xend=.end, 
                                              y=2^score, 
                                              yend=2^score), 
                                          inherit.aes=FALSE, 
                                          size=3, 
                                          col=cols[1], alpha=1, na.rm=TRUE)
              }

              p1 <- p0 +  geom_segment(data=my.dat,
                                       aes(x=.start, 
                                           xend=.end, 
                                           y=r.pred, 
                                           yend=r.pred), 
                                       col=cols[2], size=2, lty=1, alpha=1,
                                       na.rm=TRUE) + 
                  ylim(c(0,ymax)) + ylab("R")

              p2 <- p1 + scale_x_continuous(breaks = (ss2$.start + ss2$.end)/2,
                                            labels=as.character(ss2$seqnames)) + 
                  theme(panel.background=element_blank(), 
                        panel.grid.major=element_line(colour="gray90", 
                                                      size=0.3), 
                        panel.grid.minor=element_blank(), 
                        axis.text.x = element_text(angle = 45, 
                                                   hjust = 1, 
                                                   size = 10))

              if(!is.na(subclone.prev)){
                  sp <- result.dat$prev[-1]
                  dd <- subset(my.dat, p <= subclone.prev)
                  cc <- brewer.pal(8,"Set2")[1:length(sp)]
                  fill <- plyr::mapvalues(dd$p, from=sp, to=cc)
                  p3 <- p2 + geom_rect(data=dd, 
                                       aes(xmin=.start, 
                                           xmax=.end, 
                                           ymin=-Inf, 
                                           ymax=Inf), 
                                       fill=fill, 
                                       alpha=0.2, 
                                       inherit.aes = FALSE)

              }else{
                  p3 <- p2
              }

              return(p3)
          }
)

dat2gr <- function(dat) {
    dat$width <- NULL
    gr <- with(dat, GRanges(seqnames, IRanges(start, end), strand=strand))
    elementMetadata(gr) <- dat[, ! names(dat) %in% c("seqnames", 
                                                     "start", 
                                                     "end", 
                                                     "strand")]
    return(gr)
}

keep24chr <- function(gr, no.sex.chr=TRUE, is.female=FALSE) {
    hg19.chr <- seqnames(hg19.seqinfo)

    if(is.female){
        hg19.chr <- setdiff(hg19.chr, "chrY")
        seqlevels(hg19.seqinfo) <- hg19.chr
    }

    if(no.sex.chr){
        hg19.chr <- setdiff(hg19.chr, c("chrX", "chrY"))
        seqlevels(hg19.seqinfo) <- hg19.chr
    }

    pp <- gr[seqnames(gr) %in% hg19.chr ]
    seqlevels(pp) <- hg19.chr
    seqinfo(pp) <- hg19.seqinfo
    return(pp)
}

#' @title getTracks
#' @rdname getTracks
#' @aliases getTracks
#' @param p1 set 1
#' @param p2 set 2
#' @param title the title
#' @return all of the requested tracks
#' @example examples/getTracks-Ex.R
getTracks <- function(p1, p2, title="") {
    p1 <- p1+
        theme(legend.position="none",
              axis.text.x=element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin=unit(c(1,1,0,1), "cm"))

    if(title !=  ""){
        p1 <- p1 + labs(title=title) + theme(plot.title=element_text(hjust=0))
    }

    p2 <- p2 +
        theme(legend.position="none",
              axis.title.x = element_blank(),
              plot.title = element_blank(),
              plot.margin=unit(c(0,1,1,1), "cm"),
              strip.background = element_blank(),
              strip.text.x = element_blank())

    gp1 <- ggplot_gtable(ggplot_build(p1))
    gp2 <- ggplot_gtable(ggplot_build(p2))
    maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
    gp1$widths[2:3] <- maxWidth
    gp2$widths[2:3] <- maxWidth

    z <- arrangeGrob(gp1, gp2, ncol=1)
    return(z)
}

