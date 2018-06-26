test_that("VCF2VE works", {
    vcf <- SeqArray::seqExampleFileName("vcf")
    ve <- VCF2VE(vcf, out.dir = tempfile())

    ## general
    expect_equal(dim(ve), c(1348L, 90L))
    expect_equal(dim(rowData(ve)), c(1348L, 14L))
    expect_equal(dim(colData(ve)), c(90L, 0L))
    expect_s4_class(rowData(ve), "DelayedDataFrame")
    expect_s4_class(colData(ve), "DelayedDataFrame")

    ## contents
    expect_equal(names(assays(ve)),
                 c("genotype/data", "phase/data", "annotation/format/DP/data"))
    expect_equal(colnames(rowData(ve)),
                 c("ID", "ALT", "REF", "QUAL", "FILTER", "info_AA", "info_AC",
                   "info_AN", "info_DP", "info_HM2", "info_HM3", "info_OR",
                   "info_GP", "info_BN" ))
    expect_equal(colnames(colData(ve)), character(0))

    ## arguments
    ve1 <- VCF2VE(vcf, out.dir = tempfile(), info.import=c("OR", "GP")) 
    expect_equal(colnames(rowData(ve1)),
                 c("ID", "ALT", "REF", "QUAL", "FILTER", "info_OR", "info_GP")) 

    ve2 <- VCF2VE(vcf, replace=TRUE, fmt.import="xx")
    expect_equivalent(ve, ve2)
    ve3 <- VCF2VE(vcf, replace=TRUE, fmt.import=NULL)
    expect_equivalent(ve, ve3)
    ve4 <- VCF2VE(vcf, replace=TRUE, fmt.import="DP")
    expect_equivalent(ve, ve4)
    ve6 <- VCF2VE(vcf, replace=TRUE, reference="hg18")
    ## "reference" only works when vcf head does not have reference entry. 
    expect_equivalent(ve, ve6)
    
    ve7 <- VCF2VE(vcf, replace=TRUE, start=101, count=1000)  ## "start" and "count" argument works.
    expect_equal(dim(ve7), c(1000L, 90L))
    expect_equal(rownames(ve7), as.character(seq(101, length.out=1000)))
    
    ve8 <- VCF2VE(vcf, replace=TRUE, sample.info = "VE.temp/sample.info.txt")
})
