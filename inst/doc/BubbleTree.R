## ------------------------------------------------------------------------
library(BubbleTree)
load(system.file("data", "allCall.lst.RData", package="BubbleTree"))

btreeplotter <- new("BTreePlotter", max.ploidy=5, max.size=10)

pdf.filename <- paste("adjplot_all", "pdf", sep=".")
pdf(pdf.filename, width=6, height=4)

l_ply(names(allCall.lst), function(nn){
    cat("Processing sample name: ", nn, "\n")
    rbd1 <- allCall.lst[[nn]]@rbd
    rbd2 <- allCall.lst[[nn]]@rbd.adj
    arrows <- trackBTree(btreeplotter, rbd1, rbd2, min.srcSize=0.01, min.trtSize=0.01)
    btree <- drawBTree(btreeplotter, rbd1) + drawBubbles(btreeplotter, rbd2, "gray80") + arrows
    print(btree)
})
dev.off()

pdf.filename <- paste("adjplot_single", "pdf", sep=".")
pdf(pdf.filename, width=8, height=6)

l_ply(names(allCall.lst), function(nn){
    cat("Processing sample name: ", nn, "\n")
    rbd1 <- allCall.lst[[nn]]@rbd
    rbd2 <- allCall.lst[[nn]]@rbd.adj
    arrows <- trackBTree(btreeplotter, rbd1, rbd2, min.srcSize = 0.1, min.trtSize= 0.1)
    btree <- drawBTree(btreeplotter, rbd1)  + arrows
    print(btree)
})
dev.off()


## ------------------------------------------------------------------------
library(BubbleTree)
load(system.file("data", "allRBD.lst.RData", package="BubbleTree"))
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
    cat("Processing sample name: ", nn, "\n")
    rbd <- allRBD.lst[[nn]]
    btreepredictor@config$high.ploidy <- high.ploidy[nn]
    btreepredictor@config$high.purity <- high.purity[nn]
    btreepredictor <- loadRBD(btreepredictor, rbd)
    btreepredictor@config$min.segSize <- ifelse(max(btreepredictor@rbd$seg.size, na.rm=TRUE) < 0.4, 0.1, 0.4)
    
    btreepredictor <- btpredict(btreepredictor)
    
    cat(info(btreepredictor), "\n")
    return(btreepredictor)
})

names(allCall.lst) <- names(allRBD.lst)
results <- list()
for (name in names(allCall.lst)) {
    results[[name]] <- allCall.lst[[name]]@result$dist
}

xls.filename <- paste("all_calls_report", "xlsx", sep=".")
print(xls.filename)
saveXLS(results, xls.filename)


## ------------------------------------------------------------------------
library(BubbleTree)
load(system.file("data", "allCall.lst.RData", package="BubbleTree"))
load(system.file("data", "cancer.genes.minus2.rda", package="BubbleTree"))
load(system.file("data", "vol.genes.rda", package="BubbleTree"))
load(system.file("data", "gene.uni.clean.gr.rda", package="BubbleTree"))
load(system.file("data", "cyto.gr.rda", package="BubbleTree"))

# 77 common cancer genes
comm <- btcompare(vol.genes, cancer.genes.minus2)

btreeplotter <- new("BTreePlotter", branch.col="gray50")
annotator <-new("Annotate")

pdf.filename <- paste("bubbletree_preview", "pdf", sep=".")
pdf(pdf.filename, width=10, height=6)

allPreview <- llply(names(allCall.lst), function(nn){
    cat("Processing sample name: ", nn, "\n")
    cc <- allCall.lst[[nn]]
    z <- drawBTree(btreeplotter, cc@rbd.adj) + ggplot2::labs(title=sprintf("%s (%s)", nn, info(cc)))
    print(z)
    out <- cc@result$dist  %>% filter(seg.size >= 0.1 ) %>% arrange(gtools::mixedorder(as.character(seqnames)), start) # needs to be relevel
    
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


## ------------------------------------------------------------------------
library(BubbleTree)
load(system.file("data", "allCall.lst.RData", package="BubbleTree"))
load(system.file("data", "centromere.dat.rda", package="BubbleTree"))
load(system.file("data", "all.somatic.lst.RData", package="BubbleTree"))
load(system.file("data", "allHetero.lst.RData", package="BubbleTree"))
load(system.file("data", "allCNV.lst.RData", package="BubbleTree"))
load(system.file("data", "hg19.seqinfo.rda", package="BubbleTree"))

trackplotter <- new("TrackPlotter")

gr2 = centromere.dat

pdf.filename <- paste("all_tracks3", "pdf", sep=".")
pdf(pdf.filename, width=10, height=7)

plyr::l_ply(names(allCall.lst), function(nn) {
    cat("Processing sample name: ", nn, "\n")
    ymax <- ifelse(nn %in% c("lung.wgs", "lung.wes"), 9, 4.3)
    
    p1 <- xyTrack(trackplotter,
                  result.dat=allCall.lst[[nn]]@result,
                  gr2=gr2,
                  ymax=ymax) + ggplot2::labs(title=nn)
    
    p2 <- bafTrack(trackplotter,
                   result.dat=allCall.lst[[nn]]@result,
                   gr2=gr2,
                   somatic.gr=all.somatic.lst[[nn]])
    
    t1 <- getTracks(p1, p2)
    
    z1 <- heteroLociTrack(trackplotter, allCall.lst[[nn]]@result, gr2, allHetero.lst[[nn]])
    z2 <- RscoreTrack(trackplotter, allCall.lst[[nn]]@result, gr2, allCNV.lst[[nn]])
    t2 <- getTracks(z1, z2)
    
    gridExtra::grid.arrange(t1,t2, ncol=1)
    
})
dev.off()


## ------------------------------------------------------------------------
library(BubbleTree)
load(system.file("data", "allCall.lst.RData", package="BubbleTree"))

btp <- new("BTreePlotter", max.ploidy=5, max.size=10)
pairs <- list(c("lung.wes", "lung.wgs"),
              c("HCC4.Primary.Tumor", "HCC4.Recurrent.Tumor"),
              c("HCC11.Primary.Tumor", "HCC11.Recurrent.Tumor" ))

pdf("medi.comparePlot.single2.pdf", width=8,height=6)

plyr::l_ply(pairs, function(p){
    cat(p, "\n")
    rbd1 <- allCall.lst[[p[1]]]@result$dist
    rbd2 <- allCall.lst[[p[2]]]@result$dist
    
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
        z <- z + arrows + ggplot2::labs(title=sprintf("%s -> %s", p[1], p[2]))
    print(z)
})
dev.off()


## ------------------------------------------------------------------------
library(BubbleTree)
load(system.file("data", "allCall.lst.RData", package="BubbleTree"))

all.summary <- plyr::ldply(allCall.lst, function(.Object) {
    purity <- .Object@result$prev[1]
    adj <- .Object@result$ploidy.adj["adj"]
    ploidy <- (2*adj -2)/purity + 2  # when purity is low the calculation result is not reliable
    
    with(.Object@result,
         return(c(Purity=round(purity,3),
                  Prevalences=paste(round(prev,3), collapse=", "),
                  "Tumor ploidy"=round(ploidy,1))))
}) %>% plyr::rename(c(".id"="Sample"))

xls.filename <- paste("all_summary", "xlsx", sep=".")
saveXLS(list(Summary=all.summary), xls.filename)


