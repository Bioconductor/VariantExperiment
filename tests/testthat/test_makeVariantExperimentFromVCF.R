test_that("makeVariantExperimentFromVCF works", {
    vcf <- SeqArray::seqExampleFileName("vcf")
    ve <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile())
    ## general
    expect_equal(dim(ve), c(1348L, 90L))
    expect_equal(dim(rowData(ve)), c(1348L, 14L))
    expect_equal(dim(colData(ve)), c(90L, 0L))
    expect_s4_class(rowData(ve), "DelayedDataFrame")
    expect_s4_class(colData(ve), "DelayedDataFrame")
    expect_equal(rownames(rowData(ve)), rownames(ve))
    expect_equal(rownames(colData(ve)), colnames(ve))
    expect_s4_class(assay(ve,1), "DelayedArray") ## "assay,SummarizedExperiment"
                                                 ## reassign
                                                 ## dimnames[1:2] and
                                                 ## thus demotes
                                                 ## GDSArray as
                                                 ## DelayedArray.
    
    ## contents
    expect_equal(names(assays(ve)),
                 c("genotype/data", "phase/data", "annotation/format/DP/data"))
    expect_equal(colnames(rowData(ve)),
                 c("ID", "ALT", "REF", "QUAL", "FILTER", "info_AA", "info_AC",
                   "info_AN", "info_DP", "info_HM2", "info_HM3", "info_OR",
                   "info_GP", "info_BN" ))
    expect_equal(colnames(colData(ve)), character(0))

    ## arguments
    ve1 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), info.import=c("OR", "GP")) 
    expect_equal(colnames(rowData(ve1)),
                 c("ID", "ALT", "REF", "QUAL", "FILTER", "info_OR", "info_GP")) 

    ve2 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), fmt.import="xx")
    expect_equivalent(ve, ve2)
    ve3 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), fmt.import=NULL)
    expect_equivalent(ve, ve3)
    ve4 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), fmt.import="DP")
    expect_equivalent(ve, ve4)
    ve6 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(), reference="hg19")
    expect_equivalent(ve, ve6)
    ## FIXME: where is the reference info saved in VE? gds attributes? 
    
    ve7 <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(),
                                        start=101, count=1000) 
    expect_equal(dim(ve7), c(1000L, 90L))
    expect_equal(rownames(ve7), as.character(seq(101, length.out=1000)))


    sample.info <- system.file("extdata", "Example_sampleInfo.txt",
                               package="VariantExperiment")
    ve8 <- suppressWarnings(
        makeVariantExperimentFromVCF(vcf, out.dir = tempfile(),
                                     sample.info = sample.info))
    expect_equal(dim(colData(ve8)), c(90L, 1L))
    expect_equal(colnames(colData(ve8)), "family")
    ## FIXME: Warning message:
    ## In (function (node, name, val = NULL, storage = storage.mode(val),  :
    ##   Missing characters are converted to "".
})
