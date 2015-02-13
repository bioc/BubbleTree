############################################
## Functions and Dataset Descriptions for ##
## BubbleTree Package                     ##
############################################


# drawBubble Function -----------------------------------------------------

#' Draw a bubble
#'
#' @description Draw single bubble to BubbleTree plot with customized
#'  label, size and color
#'
#' @aliases addBubble
#'
#' @usage drawBubble(seg.mean, hds.median, num.mark, col, min.cex=0.3,
#'   size=1, info=NULL ,adj=0.5)
#'
#' @param seg.mean copy ratio score of the segment
#' @param hds.median median HDS score of the segment
#' @param num.mark number of the marks harbored by the segment
#' @param col color of the bubble
#' @param min.cex minimum font size
#' @param size size of the bubble to scale
#' @param info label of the bubble
#' @param adj adjusted postion of the label
#'
#' @return Plots a single bubble on the BubbleTree Plot
#'
#' @examples drawBranches(main="Demo")
#' drawBubble(0.5, 0.3, 5000, col="blue", size=2, info="PTEN", adj=-0.5)
#'
#' @export

drawBubble <- function(seg.mean, hds.median, num.mark, col, min.cex=0.3, size=1,
                       info=NULL, adj=0.5){

  col <- col
  points(2^as.numeric(seg.mean), hds.median, pch=21,
         cex=pmax(log2(num.mark/1000)*size, min.cex),
         col=col, bg=add.transparency(col, 150))

  if(!is.null(info) )  text(2^as.numeric(seg.mean), hds.median, info,
                            cex=pmax(num.mark/5000, min.cex), col=col, adj=adj)

}

# drawBranches Function ---------------------------------------------------

#' Plot branches of BubbleTree plot
#'
#' @description Plot branches of BubbleTree plot
#'
#' @usage drawBranches(xmax=3.2, main="")
#'
#' @param xmax define the upper limit of the x-axis
#' @param main title of the plot
#'
#'
#' @details The branches of BubbleTree plot stand for interger copy number change.
#'   For example, "B" and "BB" indicates LOH and copy-number neutral LOH, repectively.
#'
#' @return A plot showing branches of a BubbleTree
#'
#' @examples drawBranches(xmax=2.6)
#' drawBranches()
#'
#' @export

drawBranches <- function(xmax=3.2, main=""){
  cex=0.6
  adj=0.4
  offset=0

  p <- seq(0,1, length=11)
  p.cex <- p+0.2
  p.text <- floor(p*100)
  text.cex <- 0.5*p + 0.2
  # plot ratio and baf deviation
  # A
  plot((2-p)/2, p/(4-2*p), type='b', xlim=c(0,xmax), ylim=c(0,0.5), xlab="R score", ylab="HDS", cex=p.cex, pch="",main=main)
  text((2-p)/2, p/(4-2*p), p.text, cex=text.cex)
  text(0.5,0.5, "B", cex=cex, adj=1.9+adj, offset=-offset, col="brown")

  # homozygous deletion
  points(1-p, rep(0, length(p)), type='b', col="black", cex=p.cex, pch="")
  text(1-p, rep(0, length(p)), p.text, cex=text.cex, srt=90)
  text(0,0.02, expression(phi), cex=cex, adj=0.5, col="brown")

  # AAB,
  points(1+p/2, p/(4+2*p), type='b', col="black", cex=p.cex, pch="")
  text(1+p/2, p/(4+2*p), p.text, cex=text.cex)
  text(1.5,0.15, "ABB", cex=cex, adj=-adj-0.2, offset=offset, col="brown")
  # AAA
  points(1+p/2, 3*p/(4+2*p), type='b', col="black", cex=p.cex, pch="")
  text(1+p/2, 3*p/(4+2*p), p.text, cex=text.cex)
  text(1.5,0.5, "BBB", cex=cex, adj=-adj-0.2, offset=offset, col="brown")

  # AAAA
  points(1+p, p/(1+p), type='b', col="blue", cex=p.cex, pch="")
  text(1+p, p/(1+p), p.text, cex=text.cex, col="blue")
  text(2,0.5, "BBBB", cex=cex, adj=0-adj, offset=offset, col="blue")

  # AABB
  points(1+p, rep(0, length(p)), type='b', col="blue", cex=p.cex, pch="")
  text(1+p, rep(0, length(p)), p.text, cex=text.cex,srt=90, col="blue")
  text(2,0.02, "AABB", cex=cex, adj=0-adj, col="blue")

  # AAAB, which is similar to AAB
  points(1+p, p/(2+2*p), type='b', col="blue", pch=16,  cex=p.cex)
  text(2,0.22, "ABBB", cex=cex, adj=0-adj, col="blue", offset=offset)

  # AA
  points(rep(1, length(p)), p/2, type='b', col="black",  cex=p.cex, pch="")
  text(rep(1, length(p)), p/2, p.text, col="black",  cex=text.cex)
  text(1,.5, "BB", cex=cex, adj=1.5+adj,  offset=-offset, col="brown")

  # AB
  points(rep(1, length(p)), rep(0, length(p)), type='p', pch=16, col="red", lty=2)

  # add AAAAA, AAAAB, AAABB
  .hds <- function(x,y,p){
    p*abs(y-x)/(2*((x+y)*p  + 2*(1-p)))
  }

  .cn.ratio <- function(x,y,p){
    1-p + (x+y)*p/2
  }
  points(.cn.ratio(0,5,p), .hds(0,5,p), type='b', col="darkgreen",  cex=p.cex, pch="")
  text(.cn.ratio(0,5,p), .hds(0,5,p), p.text, col="darkgreen",  cex=text.cex)
  text(2.7,.45, "BBBBB", cex=cex, adj=1+adj,  offset=-offset, col="darkgreen")

  #ABBBB
  points(.cn.ratio(1,4,p), .hds(1,4,p), type='b', col="darkgreen",  cex=p.cex, pch="")
  text(.cn.ratio(1,4,p), .hds(1,4,p), p.text, col="darkgreen",  cex=text.cex)
  text(2.7,.32, "ABBBB", cex=cex, adj=1+adj,  offset=-offset, col="darkgreen")

  #AABBB
  points(.cn.ratio(2,3,p), .hds(2,3,p), type='b', col="darkgreen",  cex=p.cex, pch="")
  text(.cn.ratio(2,3,p), .hds(2,3,p), p.text, col="darkgreen",  cex=text.cex)
  text(2.7,.08, "AABBB", cex=cex, adj=1+adj,  offset=-offset, col="darkgreen")

  if(xmax > 3){
    # add 6x
    # AAAAAA/BBBBBB
    .n=6
    .col <- "purple"
    .x <- rep(3.2,4)
    .y <- c(0.45, 0.30, 0.14, 0.02)
    .labels <- c("BBBBBB", "ABBBBB", "AAAABB", "AAABBB")
    for(j in 0:3){
      points(.cn.ratio(j,.n-j,p), .hds(j,.n-j,p), type='b', col=.col,  cex=p.cex, pch="")
      text(.cn.ratio(j,.n-j,p), .hds(j,.n-j,p), p.text, col=.col,  cex=text.cex)
      text(.x,.y, .labels, cex=cex, adj=1+adj,  offset=-offset, col=.col)
    }

  }

}

# calc.prev function ------------------------------------------------------


#' Calculate tumor cell prevalence in a sample, an indication of sample purity
#'
#' @description A model fitting of component distributions calculated from
#' somatic copy number segments within the BubbleTree diagram.
#'
#' @aliases calc.prev
#'

#' @import mixdist
#' @importFrom geosphere distm
#' @importFrom geosphere dist2Line
#' @importFrom shape Arrows
#' @importFrom shape Arrowhead
#'
#' @usage calc.prev(rbdx,heurx=FALSE,modex=3,plotx="~/prev_model.pdf")
#'
#' @param   rbdx   a matrix or data frame generated by the function plotBubble().
#' @param   heurx  a logical value indicating if only using A/B and AA/BB
#' branches for purity calculation or all branches.
#' @param   modex  an integer value providing the expected number of modes in
#' mixture distribution.
#' @param   plotx  a character string specifying the mixture model plot file
#' name.
#'
#'
#' @details The top n (user defined) most frequent prevalence estimates are
#' seeded as means in this finite mixture model, to predict the tumor
#' (sub)clonal prevalences. This allows a user-defined expected number of
#' modes within the component distribution, though overlapping modes converge
#' to the same estimate. Then these SCNA segment frequencies are used to
#' estimate the means (+/- standard deviations) of the component distributions
#' using the R package mixdist [Macdonald et al, 1988]. The standard deviations
#' are constrained by the Poisson relation given by
#'  \eqn{\alpha=\sqrt{u_i} ,i=1,...,k}{\alpha=\sqrt(u_i) ,i=1,...,k}.  .
#'
#' @return  List object of two elements: 1) the rbdx data frame with two
#' addition columns (2^seg.mean and genotype with frequency for each segment),
#'  and 2) a data matrix of modex rows and two columns indicating the seeding
#'   modes (column 1) and estimated modes (column 2), each with the number of
#'    segments supporting each mode, separated by an underscore. The largest
#'     mode in column two is the estimated tumor purity
#'
#' @examples
#' #load sequence variants
#' data('hetero.gr', package='BubbleTree', envir = environment())
#' #load copy number variation data
#' data('cnv.gr', package='BubbleTree', envir = environment())
#' rbd<-getRBD(hetero.gr, cnv.gr) #plot BubbleTree
#' pur <- calc.prev(rbdx=rbd,heurx=FALSE,modex=3,plotx="prev_model.pdf")
#'
#'# extract the genotype (branch) and frequency for each segment
#'  pur[[1]]$ploidy_prev
# # tumor purity
#'  pur[[2]][nrow(pur[[2]]),2]
#'
#' @references
#' Macdonald PDM. and Green PEJ: User's Guide to Program MIX: An Interactive
#'  Program for Fitting Mixtures of Distributions. ICHTHUS DATA SYSTEMS 1988.
#' Macdonald, PDM (1988). \emph{Demonstration Examples for MIX 2.3}. Ichthus
#'  Data Systems, Hamilton, Ontario. 13 pp. ISBN 0-9692305-1-4.
#'
#' @export
calc.prev <- function(rbdx,heurx=FALSE,modex=3,plotx="~/prev_model.pdf") {
  # require(geosphere)
  # require(mixdist)
  .hds <- function(x,y,p){
    p*abs(y-x)/(2*((x+y)*p  + 2*(1-p)))
  }
  .cn.ratio <- function(x,y,p){
    1-p + (x+y)*p/2
  }
  .dist.pnt.line <- function(x,pnty) {
    d <- dist2Line(pnty, x)
    d
  }
  .extract.col <- function(x,colx) {
    x[,colx]
  }
  .n <- 6

  # define points (segments)
  rbdx$x <- 2^rbdx$seg.mean
  pntx <- data.frame(rbdx$x,rbdx$hds.median)

  # define lines (branches)
  p <- seq(0,1, length=41)  # prev intervals of 2.5%
  p.text <- floor(p*100)

  linex <- list(cbind((2-p)/2,p/(4-2*p)),cbind((1-p),rep(0,length(p))),
                cbind(1+p/2,p/(4+2*p)),cbind(1+p/2,3*p/(4+2*p)),
                cbind(1+p,p/(1+p)),cbind(1+p,rep(0,length(p))),cbind(1+p,p/(2+2*p)),
                cbind(rep(1, length(p)),p/2),
                cbind(.cn.ratio(0,5,p),.hds(0,5,p)),cbind(.cn.ratio(1,4,p),.hds(1,4,p)),
                cbind(.cn.ratio(2,3,p),.hds(2,3,p)),
                cbind(.cn.ratio(0,.n-0,p),.hds(0,.n-0,p)),cbind(.cn.ratio(1,.n-1,p),
                                                                .hds(1,.n-1,p)),
                cbind(.cn.ratio(2,.n-2,p),.hds(2,.n-2,p)),cbind(.cn.ratio(3,.n-3,p),
                                                                .hds(3,.n-3,p)))
  ploidy <- c("A/B","Phi","AAB/ABB","AAA/BBB","AAAA/BBBB","AABB","AAAB/ABBB",
              "AA/BB","AAAAA/BBBBB","AAAAB/ABBBB","AAABB/AABBB",
              "AAAAAA/BBBBBB","AAAAAB/ABBBBB","AABBBB/AAAABB","AAABBB")
  names(linex)<- ploidy

  # calculate distances from each segment to each branch
  tmp <- lapply(linex,.dist.pnt.line,pnty=pntx)
  distx <- matrix(unlist(lapply(tmp,.extract.col,colx=1)),nrow=nrow(tmp[[1]]),
                  ncol=length(tmp))
  x.distx <- matrix(unlist(lapply(tmp,.extract.col,colx=2)),nrow=nrow(tmp[[1]]),
                    ncol=length(tmp))
  y.distx <- matrix(unlist(lapply(tmp,.extract.col,colx=3)),nrow=nrow(tmp[[1]]),
                    ncol=length(tmp))

  distx <- as.data.frame(distx)
  x.distx <- as.data.frame(x.distx)
  y.distx <- as.data.frame(y.distx)

  dimnames(distx)[[2]] <- ploidy
  dimnames(x.distx)[[2]] <- ploidy
  dimnames(y.distx)[[2]] <- ploidy

  # for each segment, calculate ploidy and prev
  purity <- NULL
  for(i in 1:nrow(distx)) {
    branch <- names(distx)[as.logical(min(distx[i,])==distx[i,])]
    a <- NULL
    for(j in branch) {
      prev <- p.text[which(distm(c(x.distx[i,j],y.distx[i,j]),
                                            linex[[j]])==min(distm(c(x.distx[i,j],
                                                                                y.distx[i,j]),linex[[j]])))]
      a <- c(a,paste(j,"_",prev,sep=""))
    }
    a <- paste(a,collapse="//")
    purity <- rbind(purity,a)
  }

  rbdx$ploidy_prev <- purity

  if(heurx) {purity <- purity[c(grep("A/B_",purity),grep("AA/BB_",purity))]}
  y <- as.numeric(gsub("[A-Z,a-z,_,/]","",unlist(strsplit(purity,"//"))))
  prevx <- table(y)

  # plot discrete histogram and density of prev counts
  pdf(plotx)
  par(mfrow=c(2,1),oma=c(0.1,5,0.1,5))
  h <- hist(y,breaks=20,plot=FALSE)
  barplot(h$counts,names.arg=h$mids,col="lightblue",ylab="Segment counts",
          xlab="Prevalence estimates (%)",main="",cex.names=.85)
  midx <- h$counts
  names(midx) <- h$mid

  # capture the top 5 most prevalent segments for seeding model
  b <- rev(sort(midx))[1:modex]
  b <- b[order(as.numeric(names(b)))]
  priors <- as.numeric(names(b))

  # fit finite mixture dist'n w/ mle using newton-type method and em algorithms
  parx <- data.frame(prev=as.numeric(names(prevx[order(as.numeric(names(prevx)))])),
                     counts=as.numeric(prevx[order(as.numeric(names(prevx)))]))
  mix.mod1 <- mix(mixdat=parx, mixparam(mu=priors,sigma=rep(1,modex)),
                           dist="pois",emsteps=10,
                           constr=mixconstr(consigma = "POIS"),iterlim=500)
  plot.mix(mix.mod1)
  dev.off()

  k <- data.frame(prev_segcnt=paste(names(b),as.numeric(b),sep="_"),
                  prev_model=paste(round(mix.mod1$parameters[,"mu"],1),
                                   "+/-",round(mix.mod1$parameters[,"sigma"],1),sep=""))
  return(list(rbdx,k))
}

# compareBubbles function -------------------------------------------------

#'compare bubbles from two samples
#'
#' @description  TBA
#'
#' @aliases compareBubbles
#'
#'
#' @import GenomicRanges
#' @import IRanges
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom plyr dlply
#' @importFrom plyr .
#' @importFrom plyr ddply
#
#'
#' @usage compareBubbles(rbd1, rbd2, min.mark=500, min.dist=0.2, max.dist=100, main="")
#'
#' @param rbd1    RBD (R-score BAF Dataframe) from the sample 1
#' @param rbd2    RBD data.frame from the sample 2
#' @param min.mark  integer segments with mininum markers to be compared
#' @param min.dist  numeric minimum distance of the overlapped segments to be displayed
#' @param max.dist  numeric maximum distance of the overlapped segments to be displayed
#' @param main  character string for the plot title
#'
#' @details The segments (larger than min.mark) from the two samples are compared
#'  to each other.
#'
#' @return   A list of the detailed information of the overlapped segments
#'
#' @examples
#'
#' data('hcc.rbd.lst', package='BubbleTree', envir = environment())
#'
#' # show the SCNV changes between the recurrent tumor and the primary tumor
#' compareBubbles(hcc.rbd.lst$HCC11.Primary.Tumor,
#'     hcc.rbd.lst$HCC11.Recurrent.Tumor,min.dist=0.05, min.mark=2000)
#'
#' # show the similarity in the recurrent tumors between two subjects
#' # Interestingly, 17p- and 17q+ are conserved.
#' compareBubbles(hcc.rbd.lst$HCC4.Recurrent.Tumor,
#'     hcc.rbd.lst$HCC11.Recurrent.Tumor,
#'     min.dist=0.0, max.dist=0.1, min.mark=500)
#' @export


compareBubbles <- function(rbd1, rbd2,  min.mark = 500,  min.dist = 0.2,
                           max.dist = 100, main=""){

  rbd1 <- subset(rbd1, rbd1$num.mark >= min.mark)
  rbd2 <- subset(rbd2, rbd2$num.mark >= min.mark)
  pt.gr <- GRanges(rbd1$chr, IRanges(rbd1$start, rbd1$end), id=rbd1$seg.id) %>%
       annoByCytoband
  rt.gr <- GRanges(rbd2$chr, IRanges(rbd2$start, rbd2$end), id=rbd2$seg.id) %>%
       annoByCytoband
  ov <- findOverlaps(pt.gr, rt.gr, minoverlap = 1e6) %>% as.data.frame
  out <- apply(ov, 1, function(rr){
    # calculate the distance between the two
    x <- rr[1]
    y <- rr[2]
    dis <- sqrt((rbd1[x[1], "hds.median"] - rbd2[y, "hds.median"])^2 + (2^rbd1[x, "seg.mean"] - 2^rbd2[y, "seg.mean"])^2)
  })

  ov$dist <- out
  top <- ov %>% filter(dist >= min.dist & dist <= max.dist)
  if(nrow(top) == 0){
    cat("No bubble to show!\n")
    return(NULL)
  }
  ck <- chr.ck(22) # color
  drawBranches(xmax=3.2, main=main)
  rst <- apply(top, 1, function(rr){
    a <- rr[1]
    b <- rr[2]

    out <- list(pt=rbd1[a,], rt=rbd2[b,], dist = rr[3])

    x0 <- out$pt$seg.mean
    y0 <- out$pt$hds.median
    x1 <- out$rt$seg.mean
    y1 <- out$rt$hds.median

    drawBubble(x0, y0, out$pt$num.mark, ck$col[out$pt$chr],
                            info=pt.gr$cyto.band[a], size=1.5)
    drawBubble(x1, y1, out$rt$num.mark, ck$col[out$rt$chr],
                            info=rt.gr$cyto.band[b], size=1.5)



    col <- ck$col[out$pt$chr]
    Arrows(2^x0, y0, 2^x1 , y1, arr.length = 0, code = 2,
                  arr.col = "red", col=col,  lwd = 1, arr.type="curved",arr.adj=1)
    Arrowhead(2^x0, y0,  arr.length = 0.05, arr.type = "circle",
                     arr.col = col, lcol=col, arr.adj=0)
    Arrows(2^x0, y0, 2^x1 , y1, arr.length = 0.2, code = 2,
                  arr.col = add.transparency("red", 70), col="red",  lwd = 1,
                  arr.type="curved",arr.adj=1, segment=FALSE)
    return(out)
  })
  #dev.off()
  invisible(rst)
}

# annoByCytoband function -------------------------------------------------

annoByCytoband <- function(gr){
  data(cyto.gr)
  cc  <- annoByOverlap(gr=gr, gene.grs=cyto.gr, ann=cyto.gr$name)
  gr$cyto.band <- paste(sub("chr", "", seqnames(gr)), sapply(cc, function(x) {
    if(length(x)  == 1)
      return(x)
    else{
      # check whether the first and the last have the same prefix
      z <- tail(x, 1)
      if(substr(z,1,1) == substr(x[1],1,1)){
        z <- sub("p|q", "", tail(x, 1))
      }
      return(paste(x[1], z, sep='-'))
    }
  }), sep="")
  return(gr)
}

# annoByOverlap function --------------------------------------------------

annoByOverlap <- function(chr, start, end, gene.grs=NULL, gr=NULL,
                          ann=gene.grs$gene.symbol){
  if(is.null(gr)){
    gr <- GRanges(chr, IRanges(start, end))
  }
  ol <- as.data.frame(findOverlaps(gr, gene.grs))
  ol$ann <- ann[ol$subjectHits]
  rv <- dlply(ol, .(queryHits), function(df) df$ann)
  out <- list()
  out[1:length(gr)] <- ""
  out[as.numeric(names(rv))] <- rv
  return(out)
}


getCS <- function(snp.gr, cnv.gr){
  # require(GenomicRanges)
  hits = findOverlaps(snp.gr, cnv.gr)
  hits.df=as.data.frame(hits)
  snp.df=as.data.frame(snp.gr)
  cnv.df=as.data.frame(cnv.gr)
  colnames(cnv.df)=gsub("start","cnv.start",colnames(cnv.df))
  colnames(cnv.df)=gsub("end","cnv.end",colnames(cnv.df))
  cnv.df$seg.id=1:nrow(cnv.df)

  out=cbind(snp.df[hits.df$queryHits,],
            cnv.df[hits.df$subjectHits,c("num.mark","seg.mean",
                                         "seg.id","cnv.start",
                                         "cnv.end")])

  gr0=GRanges(out$seqnames, IRanges(start=out$start, end=out$end))
  hits.cy = as.data.frame(findOverlaps(gr0, cyto.gr))
  out$cyto.band=as.data.frame(cyto.gr)[hits.cy$subjectHits,"name"]
  return(out)
}

# getRBD function ---------------------------------------------------------

#' Get RBD (R-score BAF Dataframe) of the homogeneous SCNV segments
#'
#' @usage getRBD(snp.gr, cnv.gr, max.sd = 0.1)
#'
#' @param snp.gr a GRanges object containing BAF (B-allele frequency) of the germline
#' heterogenous loci
#' @param cnv.gr a GRanges object containing num.mark and seg.mean, generated
#' from the CNV call
#' @param max.sd Numeric value indicating the maximum standard
#' deviation of Homozygous Deviation Scores (HDS) within a cnv segment.
#' Segments with SD above this cutoff will be omitted.
#'
#' @details  This function merge BAF and CNV call results into one data frame.
#'  The segments with high
#' HDS variation are omitted. The RBD of the remaining "homogeneous" segments are
#'  returned.
#'
#' @return  A data frame to be called by plotBubble
#'
#' @examples
#' #load sequence variants
#' data('hetero.gr', package='BubbleTree', envir = environment())
#' #load copy number variation data
#' data('cnv.gr', package='BubbleTree', envir = environment())
#' rbd <- getRBD(snp.gr=hetero.gr, cnv.gr=cnv.gr)
#' plotBubbles(rbd, main="BubbleTree Plot")
#'
#' @export


getRBD <- function(snp.gr, cnv.gr, max.sd = 0.1){

  cs <- getCS(snp.gr, cnv.gr)
  .get.modes <- function(y, min.density=0){
    z <- density(y, from=0, to=1)
    # plot(z)
    rr <- rle(diff(z$y) > 0)
    dd <- cumsum(c(1,rr$lengths))
    ind <- seq_along(z$x) %in% dd[rr$values == FALSE]
    rv <- z$x[ ind & z$y > 1]
    return(rv)
  }

  cs <- subset(cs, !is.na(cs$num.mark) & !is.na(cs$seg.id) & !is.na(cs$freq ))
  cs.seg0 <- ddply(cs, .(seg.id), function(df){

    if(nrow(df) < 20){
      return(NULL)
    }
    mm <- .get.modes(df$freq) # needs to rename to more general names
    hds <-  NA

    hds <- abs(df$freq - 0.5)
    if(length(mm) == 1 & abs(mm[1] - 0.5) < 0.1){
      hds.median <-  abs(median(df$freq, na.rm=TRUE) - 0.5)
    }else{
      hds.median <- median(hds, na.rm=TRUE)
    }

    return(c(hds.median=hds.median, hds.sd = sd(hds, na.rm=TRUE),
                      num.mark=as.numeric(df$num.mark[1]),
                      seg.mean=as.numeric(df$seg.mean[1]),
                      chr=as.character(df$seqnames[1]), start=df$cnv.start[1],
                      end=df$cnv.end[1], cyto.band=as.character(df$cyto.band[1])))
  })

  cs.seg <- transform(cs.seg0, num.mark=as.numeric(cs.seg0$num.mark),
                      hds.median=as.numeric(cs.seg0$hds.median),
                      hds.sd = as.numeric(cs.seg0$hds.sd),
                      seg.mean=as.numeric(cs.seg0$seg.mean),
                      start=as.numeric(cs.seg0$start),
                      end=as.numeric(cs.seg0$end ))
  rbd <- subset(cs.seg, cs.seg$hds.sd < max.sd & !is.na(cs.seg$seg.id) )

  cat("Segments with high SD:\n")
  print(subset(cs.seg, cs.seg$hds.sd >= max.sd))

  return(rbd)
}

# plotBubbles function ----------------------------------------------------

#' Plot Bubbles
#'
#' @usage plotBubbles(rbd,  min.cex=0.3, show.cyto=TRUE, no.bayes=FALSE,
#'  xmax=3.2, size=1, main="BubbleTree Plot")
#'
#' @param rbd a data.frame containing tumor allele frequency
#' and segmented CNV previously generated by the getRBD() function
#' @param min.cex minimum size of bubble annotation on the plot
#' @param show.cyto   Logical; indicating if cytoband information
#' should be displayed on plot. Default is TRUE.
#' @param no.bayes Logical: contol labels
#' @param size   scaling factor that controls relative size of
#' bubbles appearing on plot
#' @param xmax   maximum value for R score plotted on x-axis
#' @param main  character string for the plot title
#'
#' @details  For each segment iteratively calculate the median and
#'  standard deviation (SD) of the homozygous deviation score (HDS) of the
#'   heterozygous-loci and filter those SCNAs with
#' high HDS variation (empirically, SD>0.2). The median of HDS,
#' the copy ratio and segment size for each homogenous SCNA are used to define
#'  the X-Y coordinates and sizes of the bubbles in the diagram.
#'
#' @return  Creates a bubbletree plot using the current graphics device and
#' returns a data.frame object containing summary information on genomic
#' regions affected by chromosome loss or gain used in the plot
#'
#' @examples
#' #load sequence variants
#' data('hetero.gr', package='BubbleTree', envir = environment())
#' #load copy number variation data
#' data('cnv.gr', package='BubbleTree', envir = environment())
#' rbd <- getRBD(snp.gr=hetero.gr, cnv.gr=cnv.gr)
#' plotBubbles(rbd)
#'
#' data('hcc.rbd.lst', package='BubbleTree', envir = environment())
#' pdf(file="hcc.bubbletree.pdf", width=8, height=6)
#' lapply(names(hcc.rbd.lst), function(sample) plotBubbles(hcc.rbd.lst[[sample]],
#'  size=2, main=sample))
#' dev.off()
#' @export
plotBubbles <- function(rbd, min.cex=0.3, show.cyto=TRUE, no.bayes=FALSE, xmax=3.2,
                        size=1,main="BubbleTree Plot"){
  plotBubbleBayes(rbd, min.cex=min.cex, show.cyto=show.cyto, no.bayes=TRUE,
                  xmax=xmax, size=size,main=main)
}

plotBubbleBayes <- function(cs.seg.ss, min.cex=0.3, show.cyto=TRUE, no.bayes=FALSE,
                            xmax=3.2, size=1,main=main){
  # using the universal color key
  ## assign color to each chromosome using function below
  uni.ck <- chr.ck(22)
  ck <- uni.ck
  ck$col <- uni.ck$col[cs.seg.ss$chr]

  # make plot - set graphic paramaters
  op <- par(mar=c(5,5,5,5+2))
  ## plot tree (background)
  drawBranches(xmax=xmax,main=main)
  ##  add data points
  points(2^as.numeric(cs.seg.ss$seg.mean), cs.seg.ss$hds.median, pch=21,
                         cex=pmax(log2(cs.seg.ss$num.mark/1000)*size, min.cex),
                         col=ck$col, bg=add.transparency(ck$col, 150))

  #need to select whether to add labels or not, add some or all
  # too noisy
  ind <- cs.seg.ss$num.mark > 0

  if(show.cyto & length(ind) > 0 ){
    labels <- ""
    if(!no.bayes)
      labels <- paste0(cs.seg.ss$cyto.band[ind],
                                       ' (', format(cs.seg.ss$mean.p, digits=2), ",
                                       ", format(cs.seg.ss$sd.p, digits=2),")")
    else
      labels <-cs.seg.ss$cyto.band[ind]
    # col
    cols <-  ck$col[ind]
    cols[cs.seg.ss$sd.p > 0.1] <- "black"
    text(2^as.numeric(cs.seg.ss$seg.mean[ind]), cs.seg.ss$hds.median[ind], labels,
                         cex=pmax(cs.seg.ss$num.mark/5000, min.cex), col=cols, adj=0.5)
  }
  add.legend(uni.ck$key, x=xmax+0.2, y=0.5, pch=16, cex=0.7)
  par(op)
  invisible()
}

##add a legend
add.legend <- function(col.key, cex=0.7, loc="topright", x=NULL, y=NULL,
                       title="Sample Key", xjust=0, yjust=1, xpd=NA,
                       pch=15, bty="n", horiz=FALSE, ...){
  if(is.null(x)) x <- loc
  legend(x,y, legend=names(col.key),
         pch=pch, cex=cex,  title=title,
         bty = bty, # no box around
         xjust= xjust, yjust=yjust, #alignment
         horiz=horiz,
         col=col.key, xpd=xpd,...)
  invisible()
}

#assign color to each chromosome
chr.ck <- function(N=22){
  chr <- paste0("chr", c(1:N, 'X', 'Y'))
  col <- rainbow(24)
  names(col) <- chr
  col <- col[1:N]
  return(list(col=col, key=col))
}

#adds transparency factor to colors defined by chr.ck function
add.transparency <- function(cols, alpha=150){
  c <- col2rgb(cols)
  new.c <- apply(c, 2, function(x){
    return(rgb(red=x[1], green=x[2], blue=x[3], alpha=alpha, maxColorValue=255))
  })
  return(new.c)
}

# Dataset cnv.gr ----------------------------------------------------------

#'@name cnv.gr
#'@aliases cnv.gr
#'@docType data
#'@title A sample dataset of tumor segmented copy number variations
#'@description A GRanges object containing segmented copy number log2 ratios
#'             and number of markers/segemnt in a tumor/normal sample pair from a patient
#'             with NSCLC
#'
#'@usage data(cnv.gr)
#'@format The format is: Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
#'
#'@details Metadata columns include
#'\itemize{
#'  \item \strong{num.mark} number of markers within segment
#'  \item \strong{seg.mean} mean log2 copy number ratio within a segment
#'  }
#'@examples data(cnv.gr)
#'@keywords datasets
NULL

# Dataset hetero.gr -------------------------------------------------------

#'@name hetero.gr
#'@aliases hetero.gr
#'@docType data
#'@title BAF of the germline heterozygous loci in GRanges format
#'@description A sample data of B-allele frquencies of the germline heterozygous loci
#'
#'@usage data(hetero.gr)
#'@format  The format is: Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
#'
#'@details Metadata columns include
#'\itemize{
#'  \item \strong{freq} B allele frequency
#'  }
#'@examples data(hetero.gr)
#'@keywords datasets
NULL


# Dataset hcc.rbd.lst -----------------------------------------------------

#'@name hcc.rbd.lst
#'@aliases hcc.rbd.lst
#'@docType data
#'@title Rbd data from the HCC samples
#'@description A list of RBD (R-BAF Dataframe) of the four HCC samples:\cr
#'   HCC4.Primary.Tumor \cr
#'   HCC4.Recurrent.Tumor\cr
#'   HCC11.Primary.Tumor\cr
#'   HCC11.Recurrent.Tumor\cr
#'
#'@usage data(hcc.rbd.lst)
#'@format
#'  The format is: A list of 4 data frames
#'
#'@details  These data frames are produced by the getRBD() function and include columns for:
#'\itemize{
#'   \item \strong{seg.id} copy number segment identifier
#'   \item \strong{hds.median}  median homozygosity deviation score withing the segment
#'   \item \strong{hds.sd} standard deviation of hds within the segment
#'   \item \strong{num.mark} number of markers within the segment
#'   \item \strong{seg.mean} mean copy number of markers within segment
#'   \item \strong{chr} segment chromosome
#'   \item \strong{start} segment start coordinate
#'   \item \strong{end} segment end coordinate
#'   \item \strong{cytoband} segment cytoband
#'  }
#'@examples data(hcc.rbd.lst)
#'@keywords datasets
NULL


