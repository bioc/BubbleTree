utils::globalVariables(c("seg.size", "filename", "prev", "mixedorder", 
                         "allCall.lst", "allRBD.lst", "vol.genes", 
                         "cancer.genes.minus2", "gene.uni.clean.gr", "cyto.gr",
                         "all.somatic.lst", "allHetero.lst", "allCNV.lst"))

# adjust_plot Function -----------------------------------------------------

#' Draw a BubblePlot from a list of pre-processed segments
#'
#' @description Draw a Bubble plot using pre-processed segment data that 
#' contains the following:
#'
#' seg.mean - copy ratio score of the segment
#' hds.median - median HDS score of the segment
#' num.mark - number of the marks harbored by the segment
#' col - color of the bubble
#' min.cex - minimum font size
#' size - size of the bubble to scale
#' info - label of the bubble
#' adj - adjusted postion of the label
#'
#' @aliases drawBubblePlot
#'
#' @usage adjust_plot()
#'
#' @return Plots a full BubblePlot using pre-processed segment data
#'
#' @export
adjust_plot <- function() {

    load('data/allCall.lst.RData')

    btreeplotter <- new("BTreePlotter", max.ploidy=5, max.size=10)

    pdf.filename <- paste("adjplot_all", "pdf", sep=".")
    pdf(pdf.filename, width=6, height=4)

    l_ply(names(allCall.lst), function(nn){
        cat(nn, "\n")
        rbd1 <- allCall.lst[[nn]]@rbd
        rbd2 <- allCall.lst[[nn]]@rbd.adj
        arrows <- trackBTree(btreeplotter, 
                             rbd1, 
                             rbd2, 
                             min.srcSize=0.01, 
                             min.trtSize=0.01)
        btree <- drawBTree(btreeplotter, rbd1) + 
            drawBubbles(btreeplotter, rbd2, "gray80") + 
            arrows
        print(btree)
    })
    dev.off()

    pdf.filename <- paste("adjplot_single", "pdf", sep=".")
    pdf(pdf.filename, width=8, height=6)

    l_ply(names(allCall.lst), function(nn){
        cat(nn, "\n")
        rbd1 <- allCall.lst[[nn]]@rbd
        rbd2 <- allCall.lst[[nn]]@rbd.adj
        arrows <- trackBTree(btreeplotter, 
                             rbd1, 
                             rbd2, 
                             min.srcSize = 0.1, 
                             min.trtSize= 0.1)
        btree <- drawBTree(btreeplotter, rbd1)  + arrows
        print(btree)
    })
    dev.off()
}

# call_all Function -----------------------------------------------------

#' Draw a BubblePlot from a list of pre-processed segments
#'
#' @description Draw a Bubble plot using pre-processed segment data 
#' that contains the following:
#'
#' seg.mean - copy ratio score of the segment
#' hds.median - median HDS score of the segment
#' num.mark - number of the marks harbored by the segment
#' col - color of the bubble
#' min.cex - minimum font size
#' size - size of the bubble to scale
#' info - label of the bubble
#' adj - adjusted postion of the label
#'
#' @aliases drawAllPlots
#'
#' @usage call_all()
#'
#' @return Plots a full BubblePlot using pre-processed segment data
#'
#' @export
call_all <- function() {

    load("data/allRBD.lst.RData")

    btreepredictor <- new("BTreePredictor")
    btreepredictor@config$cutree.h <- 0.15

    high.ploidy <- rep(TRUE, length(allRBD.lst))
    high.purity <- rep(TRUE, length(allRBD.lst))

    high.ploidy[c("sam6",
                  "ovary.wgs",
                  "ovary.wes",
                  "TCGA-06-0145-01A-01W-0224-08",
                  "TCGA-13-1500-01A-01D-0472-01",
                  "TCGA-AO-A0JJ-01A-11W-A071-09")] <- FALSE

    high.purity[c("sam6", "ovary.wgs", "ovary.wes")] <- FALSE


    allCall.lst <- plyr::llply(names(allRBD.lst), function(nn) {
        cat("Sample Name: ", nn, "\n")
        rbd <- allRBD.lst[[nn]]
        btreepredictor@config$high.ploidy <- high.ploidy[nn]
        btreepredictor@config$high.purity <- high.purity[nn]
        btreepredictor <- loadRBD(btreepredictor, rbd)
        btreepredictor@config$min.segSize <- 
            ifelse(max(btreepredictor@rbd$seg.size, na.rm=TRUE) < 0.4, 0.1, 0.4)

        btreepredictor <- btpredict(btreepredictor)

        cat(info(btreepredictor), "\n")
        return(btreepredictor)
    })

    names(allCall.lst) <- names(allRBD.lst)
    print(names(allCall.lst))
    results <- list()
    for (name in names(allCall.lst)) {
        results[[name]] <- allCall.lst[[name]]@result$dist
    }

    xls.filename <- paste("reports/", 
                          paste("all_calls_report", "xlsx", sep="."), sep="")
    print(xls.filename)
    saveXLS(results, xls.filename)
}

# call_default Function -----------------------------------------------------

#' A testing function that draws a BubblePlot from a list of pre-processed 
#' segments
#'
#' @description A testing function that draws a Bubble plot using pre-processed 
#' segment data that contains the following:
#'
#' seg.mean - copy ratio score of the segment
#' hds.median - median HDS score of the segment
#' num.mark - number of the marks harbored by the segment
#' col - color of the bubble
#' min.cex - minimum font size
#' size - size of the bubble to scale
#' info - label of the bubble
#' adj - adjusted postion of the label
#'
#' @aliases drawTestPlots
#'
#' @usage call_default()
#'
#' @return Plots a full BubblePlot using pre-processed segment testing data
#'
#' @export
call_default <- function() {

    load("data/allRBD.lst.RData")

    btreepredictor <- new("BTreePredictor")
    btreepredictor@config$cutree.h <- 0.15

    min.segSize <- rep(0.5, length(allRBD.lst))
    high.ploidy <- rep(TRUE, length(allRBD.lst))
    high.purity <- rep(TRUE, length(allRBD.lst))
    
    names(min.segSize) <- 
        names(high.ploidy) <- 
        names(high.purity) <- 
        names(allRBD.lst)
    
    min.segSize[c("lung.wgs", "ovary.wgs")] <- 0.1

    defaultCall.lst <- plyr::llply(names(allRBD.lst), function(nn) {
       cat(nn, "\n")

        rbd <- allRBD.lst[[nn]]
        btreepredictor@config$min.segSize <- min.segSize[nn]
        btreepredictor@config$high.ploidy <- high.ploidy[nn]
        btreepredictor@config$high.purity <- high.purity[nn]
        btreepredictor <- loadRBD(btreepredictor, rbd)
        
        btreepredictor@config$min.segSize <- 
            ifelse(max(btreepredictor@rbd$seg.size, na.rm=TRUE) < 0.4, 0.1, 0.4)
        
        btreepredictor <- btpredict(btreepredictor)
        cat(info(btreepredictor), "\n")
        print(btreepredictor)
        return(btreepredictor)
    })

    # default preview
    btreeplotter <-new("BTreePlotter", branch.col="gray50")

    pdf.filename <- paste("default_calls", "pdf", sep=".")
    pdf(pdf.filename, width=10, height=6)

    defaultPreview <- llply(names(defaultCall.lst), function(nn){
        cat(nn, "\n")
        cc <- defaultCall.lst[[nn]]
        cat("CC: ", cc@rbd.adj, "\n")
        z <- drawBTree(btreeplotter, cc@rbd.adj) + 
            labs(title=sprintf("%s (%s)", nn, info(cc)))
        print(z)
        cc@result$dist %>% 
            filter(seg.size >= 0.1 ) %>% 
            arrange(mixedorder(as.character(seqnames)), start)
    })
    dev.off()

    names(defaultPreview) <-  names(defaultCall.lst)

    names(allCall.lst) <- names(allRBD.lst)

    results <- list()
    for (name in names(allCall.lst)) {
        results[[name]] <- allCall.lst[[name]]@result$dist
    }

    xls.filename <- paste(filename, "xlsx", sep=".")
    saveXLS(results, xls.filename)

    return(length(results))
}

# plot_bubbletree_preview Function -------------------------------------------

#' A testing function that draws a BubblePlot from a list of pre-processed 
#' segments with overlays of different important cancer genes.
#'
#' @description A testing function that draws a BubblePlot from a list of 
#' pre-processed segments with overlays of different important cancer genes.
#'
#' seg.mean - copy ratio score of the segment
#' hds.median - median HDS score of the segment
#' num.mark - number of the marks harbored by the segment
#' col - color of the bubble
#' min.cex - minimum font size
#' size - size of the bubble to scale
#' info - label of the bubble
#' adj - adjusted postion of the label
#'
#' @aliases drawBubblePlotPreview
#'
#' @usage plot_bubbletree_preview()
#'
#' @return A BubblePlot with called CNVs with cancer genes overlayed
#'
#' @export
plot_bubbletree_preview <- function() {
    load("data/allCall.lst.RData")
    load("data/cancer.genes.minus2.rda")
    load("data/vol.genes.rda")
    load("data/gene.uni.clean.gr.rda")
    load("data/cyto.gr.rda")

    # 77 common cancer genes
    comm <- btcompare(vol.genes, cancer.genes.minus2)

    btreeplotter <- new("BTreePlotter", branch.col="gray50")
    annotator <-new("Annotate")

    pdf.filename <- paste("bubbletree_preview", "pdf", sep=".")
    pdf(pdf.filename, width=10, height=6)

    allPreview <- llply(names(allCall.lst), function(nn){
        cat(nn, "\n")
        cc <- allCall.lst[[nn]]
        z <- drawBTree(btreeplotter, cc@rbd.adj) + 
            labs(title=sprintf("%s (%s)", nn, info(cc)))
        
        print(z)
        out <- cc@result$dist  %>% 
            filter(seg.size >= 0.1 ) %>% 
            arrange(mixedorder(as.character(seqnames)), start) 

        ann <- with(out, {
            annoByGenesAndCyto(annotator,
                               as.character(out$seqnames),
                               as.numeric(out$start),
                               as.numeric(out$end),
                               comm$comm,
                               gene.uni.clean.gr=gene.uni.clean.gr,
                               cyto.gr=cyto.gr)
        })

        out$cyto <- ann$cyto
        out$genes <- ann$ann
        return(out)
    })
    dev.off()
}

# plot_tracks Function -----------------------------------------------------

#' A testing function that draws a BubblePlot from a list of pre-processed 
#' segments with overlays of different types of annotation data.
#'
#' @description A testing function that draws a BubblePlot from a list of 
#' pre-processed segments with overlays of different types of annotation data.
#'
#' seg.mean - copy ratio score of the segment
#' hds.median - median HDS score of the segment
#' num.mark - number of the marks harbored by the segment
#' col - color of the bubble
#' min.cex - minimum font size
#' size - size of the bubble to scale
#' info - label of the bubble
#' adj - adjusted postion of the label
#'
#' @aliases plotAllTracks
#'
#' @usage plot_tracks()
#'
#' @return A BubblePlot with called CNVs with annotation overlayed
#'
#' @export
plot_tracks <- function() {

    load("data/allCall.lst.RData")
    load("data/centromere.dat.rda")
    load("data/all.somatic.lst.RData")
    load("data/allHetero.lst.RData")
    load("data/allCNV.lst.RData")
    load("data/hg19.seqinfo.rda")

    trackplotter <- new("TrackPlotter")

    gr2 = centromere.dat

    pdf.filename <- paste("all_tracks3", "pdf", sep=".")
    pdf(pdf.filename, width=10, height=7)

    plyr::l_ply(names(allCall.lst), function(nn) {
        cat(nn, "\n")
        ymax <- ifelse(nn %in% c("lung.wgs", "lung.wes"), 9, 4.3)

        p1 <- xyTrack(trackplotter,
                      result.dat=allCall.lst[[nn]]@result,
                      gr2=gr2,
                      ymax=ymax) + labs(title=nn)

        p2 <- bafTrack(trackplotter,
                       result.dat=allCall.lst[[nn]]@result,
                       gr2=gr2,
                       somatic.gr=all.somatic.lst[[nn]])

        t1 <- getTracks(p1, p2)

        z1 <- heteroLociTrack(trackplotter, 
                              allCall.lst[[nn]]@result, 
                              gr2, 
                              allHetero.lst[[nn]])
        
        z2 <- RscoreTrack(trackplotter, 
                          allCall.lst[[nn]]@result, 
                          gr2, 
                          allCNV.lst[[nn]])
        
        t2 <- getTracks(z1, z2)

        gridExtra::grid.arrange(t1,t2, ncol=1)

    })
    dev.off()
}

# show_summary Function -----------------------------------------------------

#' Create a summary Excel file
#'
#' @description Create a summary Excel file
#'
#' @aliases showSummaryData
#'
#' @usage show_summary()
#'
#' @return Excel file with summary data
#'
#' @export
show_summary <- function() {

    load("data/allCall.lst.RData")

    all.summary <- plyr::ldply(allCall.lst, function(.Object) {
        purity <- .Object@result$prev[1]
        adj <- .Object@result$ploidy.adj["adj"]
        # when purity is low the calculation result is not reliable
        ploidy <- (2*adj -2)/purity + 2  

        with(.Object@result,
             return(c(Purity=round(purity,3),
                      Prevalences=paste(round(prev,3), collapse=", "),
                      "Tumor ploidy"=round(ploidy,1))))
    }) %>% plyr::rename(c(".id"="Sample"))

    xls.filename <- paste("all_summary", "xlsx", sep=".")
    saveXLS(list(Summary=all.summary), xls.filename)
}

# medi_compare Function -----------------------------------------------------

#' Compare BubblePlots from 2 different sets of cancer data
#'
#' @description Compare BubblePlots from 2 different sets of cancer data
#'
#' @aliases compareSamples
#'
#' @usage medi_compare()
#'
#' @return Comparison BubblePlots
#'
#' @export
medi_compare <- function() {

    load("data/allCall.lst.RData")

    btp <- new("BTreePlotter", max.ploidy=5, max.size=10)
    pairs <- list(c("lung.wes", "lung.wgs"),
                  c("HCC4.Primary.Tumor", "HCC4.Recurrent.Tumor"),
                  c("HCC11.Primary.Tumor", "HCC11.Recurrent.Tumor" ))

    #pdf("reports/medi.comparePlot.single2.pdf", width=8,height=6)
    pdf("medi.comparePlot.single2.pdf", width=8,height=6)

    plyr::l_ply(pairs, function(p){
        cat(p, "\n")
        rbd1 <- allCall.lst[[p[1]]]@result$dist
        rbd2 <- allCall.lst[[p[2]]]@result$dist

        #srcSize <- ifelse(p[1] == "lung.wgs", 0.1, 1)
        srcSize <- 0.5
        trtSize <- ifelse(p[1] == "lung.wes", 0.1, 1)
        minOver <- ifelse(p[1] == "lung.wes", 5e6, 1e7)
        arrows <- trackBTree(btp,
                             rbd1,
                             rbd2,
                             min.srcSize=srcSize,
                             min.trtSize=trtSize,
                             min.overlap=minOver)

        z <- drawBTree(btp, rbd1)
        if(!is.null(arrows))
            z <- z + arrows + labs(title=sprintf("%s -> %s", p[1], p[2]))
        print(z)
    })
    dev.off()
}
