
# load sample files
load(system.file("data", "cnv.gr.rda", package="BubbleTree"))
load(system.file("data", "snp.gr.rda", package="BubbleTree"))

# load annotations
load(system.file("data", "centromere.dat.rda", package="BubbleTree"))
load(system.file("data", "cyto.gr.rda", package="BubbleTree"))
load(system.file("data", "cancer.genes.minus2.rda", package="BubbleTree"))
load(system.file("data", "vol.genes.rda", package="BubbleTree"))
load(system.file("data", "gene.uni.clean.gr.rda", package="BubbleTree"))


# initialize RBD object
r <- new("RBD", unimodal.kurtosis=-0.1)

# create new RBD object with GenomicRanges objects for SNPs and CNVs
rbd <- makeRBD(r, snp.gr, cnv.gr)
head(rbd)

# create a new prediction
btreepredictor <- new("BTreePredictor", rbd=rbd, max.ploidy=6, prev.grid=seq(0.2,1, by=0.01))
pred <- btpredict(btreepredictor)

# create rbd plot
btreeplotter <- new("BTreePlotter", max.ploidy=5, max.size=10)
btree <- drawBTree(btreeplotter, pred@rbd)
print(btree)

# create rbd.adj plot
btreeplotter <- new("BTreePlotter", branch.col="gray50")
btree <- drawBTree(btreeplotter, pred@rbd.adj)
print(btree)

# create a combined plot with rbd and rbd.adj that shows the arrows indicating change
# THIS IS VERY MESSY WITH CURRENT DATA from Dong
btreeplotter <- new("BTreePlotter", max.ploidy=5, max.size=10)
arrows <- trackBTree(btreeplotter,
                     pred@rbd,
                     pred@rbd.adj,
                     min.srcSize=0.01, 
                     min.trtSize=0.01)

btree <- drawBTree(btreeplotter, pred@rbd) + arrows 
print(btree)


# create a plot with overlays of significant genes
btreeplotter <- new("BTreePlotter", branch.col="gray50")
annotator <- new("Annotate")

comm <- btcompare(vol.genes, cancer.genes.minus2)

sample.name <- "22_cnv_snv"

btree <- drawBTree(btreeplotter, pred@rbd.adj) + 
    ggplot2::labs(title=sprintf("%s (%s)", sample.name, info(pred)))

out <- pred@result$dist  %>% 
    filter(seg.size >= 0.1 ) %>% 
    arrange(gtools::mixedorder(as.character(seqnames)), start)

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

btree <- btree + drawFeatures(btreeplotter, out)
print(btree)


# print out purity and ploidy values
info <- info(pred)
cat("\nPurity/Ploidy: ", info, "\n")
