
load(system.file("data", "cnv.gr.rda", package="BubbleTree"))
load(system.file("data", "snp.gr.rda", package="BubbleTree"))

# initialize RBD object
r <- new("RBD", unimodal.kurtosis=-0.1)

# create new RBD object with GenomicRanges objects for SNPs and CNVs
rbd <- makeRBD(r, snp.gr, cnv.gr)

# plot a new BubbleTree plot
btreeplotter <- new("BTreePlotter", max.ploidy=5, max.size=10)
btree <- drawBTree(btreeplotter, rbd)
print(btree)


