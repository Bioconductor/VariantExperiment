gdsfile <- system.file("extdata/test.gds", package = "VariantExperiment")

test_that("grange generation works", {
    .granges_generalgds <- VariantExperiment:::.granges_generalgds
    gr <- .granges_generalgds(gdsfile, feature.num = 20)
    expect_s4_class(gr, "GRanges")
    expect_s4_class(ranges(gr), "IRanges")
})

test_that("rowRanges mcols works", {
    .rowRanges_generalgds <- VariantExperiment:::.rowRanges_generalgds
    rr <- .rowRanges_generalgds(gdsfile, "feature.id", NULL, TRUE)
    expect_s4_class(rr, "GRanges")
    expect_s4_class(mcols(rr), "DelayedDataFrame")
    expect_equal(names(mcols(rr)), c("chromosome", "position"))

    rr <- .rowRanges_generalgds(gdsfile, "feature.id", NULL, FALSE)
    expect_s4_class(rr, "GRanges")
    expect_equal(names(mcols(rr)), c("chromosome", "position"))
    expect_s4_class(mcols(rr), "DataFrame")
})

test_that("sample related nodes works", {
    .colData_generalgds <- VariantExperiment:::.colData_generalgds
    res <- .colData_generalgds(gdsfile, "sample.id", NULL, TRUE)
    expect_s4_class(res, "DelayedDataFrame")    
    expect_equal(dim(res), c(10L, 1L))
    expect_s4_class(res[[1]], "GDSArray")

    res <- .colData_generalgds(gdsfile, "sample.id", NULL, FALSE)
    expect_s4_class(res, "DataFrame")    
    expect_equal(dim(res), c(10L, 1L))
    expect_true(is.character(res[[1]]))
})

test_that("conversion function works", {
    ve <- makeVariantExperimentFromGDS(gdsfile, ftnode = "feature.id",
                                       smpnode = "sample.id")
    expect_s4_class(ve, "VariantExperiment")
    expect_equal(dim(ve), c(20L, 10L))
    expect_s4_class(assay(ve), "DelayedMatrix")
})
