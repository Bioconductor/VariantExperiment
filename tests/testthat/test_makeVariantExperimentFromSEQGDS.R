file <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")

test_that("grange generation works", {
    .granges_seqgds <- VariantExperiment:::.granges_seqgds
    gr <- .granges_seqgds(file)
    expect_s4_class(gr, "GRanges")
    expect_s4_class(ranges(gr), "IRanges")
})

test_that("varnode inmem to DF", {
    .varnode_seqgds_inmem <- VariantExperiment:::.varnode_seqgds_inmem
    res <- .varnode_seqgds_inmem(file, "allele")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), c("REF", "ALT")) ## customized for SEQ_ARRAY
})

test_that("info input works", {
    .infoColumns_seqgds <- VariantExperiment:::.infoColumns_seqgds
    res <- .infoColumns_seqgds(file, "variant.id", NULL, TRUE)
    expect_s4_class(res, "DelayedDataFrame")
    expect_equal(dim(res), c(1348L, 8L))
    expect_true(all(unlist(lapply(res, function(x) is(x, "GDSArray")))))

    res <- .infoColumns_seqgds(file, "variant.id", c("AC", "AN", "DP"), FALSE)
    expect_s4_class(res, "DataFrame")
    expect_equal(dim(res), c(1348L, 3L))
    expect_equal(colnames(res), paste0("info.", c("AC", "AN", "DP")))
    
    res <- .infoColumns_seqgds(file, "variant.id", character(0), TRUE)
    expect_equal(dim(res), c(1348L, 0L))
    res <- .infoColumns_seqgds(file, "variant.id", character(0), FALSE)
    expect_equal(dim(res), c(1348L, 0L))

    expect_error(suppressWarning(.infoColumns_seqgds(file, "variant.id", c("random"), TRUE)))
    expect_error(suppressWarning(.infoColumns_seqgds(file, "variant.id", c("random"), FALSE)))
})

test_that("rowRanges works", {
    .rowRanges_seqgds <- VariantExperiment:::.rowRanges_seqgds
    rr <- .rowRanges_seqgds(file, "variant.id", NULL, TRUE)
    expect_s4_class(rr, "GRanges")
    expect_s4_class(mcols(rr), "DelayedDataFrame")
    expect_equal(names(mcols(rr)), c(paste0("annotation.", c("id", "qual", "filter")), "REF", "ALT"))

    rr <- .rowRanges_seqgds(file, "variant.id", NULL, FALSE)
    expect_s4_class(mcols(rr), "DataFrame")
})

test_that("sample related nodes works", {
    .colData_seqgds <- VariantExperiment:::.colData_seqgds
    res <- .colData_seqgds(file, "sample.id", NULL, TRUE)
    expect_s4_class(res, "DelayedDataFrame")    
    expect_equal(dim(res), c(90L, 1L))
    expect_s4_class(res[[1]], "GDSArray")

    res <- .colData_seqgds(file, "sample.id", NULL, FALSE)
    expect_s4_class(res, "DataFrame")    
    expect_equal(dim(res), c(90L, 1L))
    expect_true(is.character(res[[1]]))
})

test_that("conversion function works", {
    ve <- makeVariantExperimentFromSEQGDS(file)
    expect_s4_class(ve, "VariantExperiment")
    expect_equal(dim(ve), c(1348L, 90L))
    expect_s4_class(assay(ve, 1), "DelayedArray")
    expect_s4_class(assay(ve, 2), "DelayedMatrix")
    expect_s4_class(assay(ve, 3), "DelayedMatrix")
    expect_s4_class(rowData(ve), "DelayedDataFrame")
    expect_s4_class(colData(ve), "DelayedDataFrame")
})

