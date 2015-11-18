load(system.file("inst/extdata", "cnv.test.gr", package="BubbleTree"))
load(system.file("inst/extdata", "snv.test.gr", package="BubbleTree"))

# initialize RBD object
r <- new("RBD", unimodal.kurtosis=-0.1)

# create new RBD object with GenomicRanges objects for SNPs and CNVs
rbd <- makeRBD(r, snp.gr, cnv.gr)
