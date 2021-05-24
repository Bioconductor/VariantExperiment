file <- SNPRelate::snpgdsExampleFileName()

test_that("grange generation works", {
    .granges_snpgds <- VariantExperiment:::.granges_snpgds
    gr <- .granges_snpgds(file)
    expect_s4_class(gr, "GRanges")
    expect_s4_class(ranges(gr), "IRanges")
})

test_that("alleles for snpgds works", {
    .varnode_snpgds_inmem <- VariantExperiment:::.varnode_snpgds_inmem
    res <- .varnode_snpgds_inmem(file, "snp.allele")
    expect_s4_class(res$snp.allele1, "DNAStringSet")
    expect_s4_class(res$snp.allele2, "DNAStringSet")
    expect_equal(dim(res), c(9088L, 2L))
})

test_that("rowRanges works", {
    .rowRanges_snpgds <- VariantExperiment:::.rowRanges_snpgds
    rr <- .rowRanges_snpgds(file, "snp.id", NULL, TRUE)
    expect_s4_class(rr, "GRanges")
    expect_s4_class(mcols(rr), "DelayedDataFrame")
    expect_equal(names(mcols(rr)), paste0("snp.", c("rs.id", "allele1", "allele2")))
    
    rr <- .rowRanges_snpgds(file, "snp.id", NULL, FALSE)
    expect_s4_class(mcols(rr), "DataFrame")
})

test_that("sample related nodes works", {
    .colData_snpgds <- VariantExperiment:::.colData_snpgds
    res <- .colData_snpgds(file, "sample.id", NULL, TRUE)
    expect_s4_class(res, "DelayedDataFrame")    
    expect_equal(dim(res), c(279L, 5L))
    expect_s4_class(res[[1]], "GDSArray")

    res <- .colData_snpgds(file, "sample.id", NULL, FALSE)
    expect_s4_class(res, "DataFrame")    
    expect_equal(dim(res), c(279L, 5L))
    expect_true(is.character(res[[1]]))
})

test_that("conversion function works", {
    ve <- makeVariantExperimentFromSNPGDS(file)
    expect_s4_class(ve, "VariantExperiment")
    expect_equal(dim(ve), c(9088L, 279L))
    expect_s4_class(assay(ve), "DelayedArray")
    expect_s4_class(rowData(ve), "DelayedDataFrame")
    expect_s4_class(colData(ve), "DelayedDataFrame")
})
