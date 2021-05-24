vcf <- SeqArray::seqExampleFileName("vcf")

test_that("makeVariantExperimentFromVCF works", {
    outdir <- tempfile()
    ve <- makeVariantExperimentFromVCF(vcf, out.dir = outdir)
    ## general
    expect_equal(dim(ve), c(1348L, 90L))
    expect_equal(dim(rowData(ve)), c(1348L, 13L))
    ## info.AA has inconsistent dim: 1328*13, so not included
    expect_equal(dim(colData(ve)), c(90L, 0L))
    expect_s4_class(rowData(ve), "DelayedDataFrame")
    expect_s4_class(colData(ve), "DelayedDataFrame")
    expect_equal(rownames(rowData(ve)), rownames(ve))
    expect_equal(rownames(colData(ve)), colnames(ve))
    expect_s4_class(assay(ve,1), "DelayedArray") ## permuted
    expect_s4_class(assay(ve,2), "DelayedMatrix")
    expect_s4_class(assay(ve,3), "DelayedMatrix")   
    
    ## ## contents
    expect_equal(names(assays(ve)),
                 c("genotype/data", "phase/data", "annotation/format/DP/data"))
    expect_equal(colnames(rowData(ve)),
                 c(paste0("annotation.", c("id", "qual", "filter")),
                          "REF", "ALT", "info.AC", "info.AN",
                          "info.DP", "info.HM2", "info.HM3",
                          "info.OR", "info.GP", "info.BN" ))
    expect_equal(colnames(colData(ve)), character(0))
    
    ## ## arguments
    ve1 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), info.import=c("OR", "GP")) 
    expect_equal(colnames(rowData(ve1)),
                 c(paste0("annotation.", c("id", "qual", "filter")), "REF",
                   "ALT", "info.OR", "info.GP")) 
    ve2 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), fmt.import="xx")
    expect_false(all(grepl("annotation/format", assayNames(ve2))))
    ve3 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), fmt.import=NULL)
    expect_equivalent(ve, ve3)
    ve4 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), fmt.import="DP")
    expect_equivalent(ve, ve4)
    ve5 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), reference="hg19")
    expect_equivalent(ve, ve5)
    ## FIXME: where is the reference info saved in VE? gds attributes? 
    
    ve6 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(),
                                        start=101, count=1000) 
    expect_equal(dim(ve6), c(1000L, 90L))
    expect_equal(rownames(ve6), as.character(seq(101, length.out=1000)))

    sample.info <- system.file("extdata", "Example_sampleInfo.txt",
                               package="VariantExperiment")
    ve7 <- suppressWarnings(
        makeVariantExperimentFromVCF(vcf, out.dir = tempfile(),
                                     sample.info = sample.info))
    expect_equal(dim(colData(ve7)), c(90L, 1L))
    expect_equal(colnames(colData(ve7)), "family")
    ## ## FIXME: Warning message:
    ## ## In (function (node, name, val = NULL, storage = storage.mode(val),  :
    ## ##   Missing characters are converted to "".
})
