
load(system.file("data", "allCall.lst.RData", package="BubbleTree"))
load(system.file("data", "cancer.genes.minus2.rda", package="BubbleTree"))
load(system.file("data", "vol.genes.rda", package="BubbleTree"))
load(system.file("data", "gene.uni.clean.gr.rda", package="BubbleTree"))
load(system.file("data", "cyto.gr.rda", package="BubbleTree"))

comm <- btcompare(vol.genes, cancer.genes.minus2)
btreeplotter <- new("BTreePlotter", branch.col="gray50")
annotator <-new("Annotate")
nn <- "sam2"
cc <- allCall.lst[[nn]]
z <- drawBTree(btreeplotter, cc@rbd.adj) + 
    ggplot2::labs(title=sprintf("%s (%s)", nn, info(cc)))
out <- cc@result$dist  %>% 
    filter(seg.size >= 0.1 ) %>% 
    arrange(gtools::mixedorder(as.character(seqnames)), start)

ann <- annoByGenesAndCyto(annotator,
                   as.character(out$seqnames),
                   as.numeric(out$start),
                   as.numeric(out$end),
                   comm$comm,
                   gene.uni.clean.gr=gene.uni.clean.gr,
                   cyto.gr=cyto.gr)
