test_that("initiate seq gds file works", {
    .initiate_seqgds <- VariantExperiment:::.initiate_seqgds
    file <- SeqArray::seqExampleFileName("gds")
    se <- makeVariantExperimentFromGDS(file)
    gds_path <- tempfile(fileext=".gds")
    .initiate_seqgds(se, gds_path, compress="LZMA_RA")
    f <- gdsfmt::openfn.gds(gds_path)
    on.exit(gdsfmt::closefn.gds(f))
    expect_true(validObject(f))
    expect_true(is(f, "gds.class"))
    expect_true(all(c("sample.id", "variant.id", "chromosome", "position") %in% gdsfmt::ls.gdsn(f)))
})

test_that("initiate snp gds file works", {
    .initiate_snpgds <- VariantExperiment:::.initiate_snpgds
    file <- SNPRelate::snpgdsExampleFileName()
    se <- makeVariantExperimentFromGDS(file)
    ## FIXME: warning
    ##   1: In methods:::.selectDotsMethod(classes, .MTable, .AllMTable) :
    ## multiple direct matches: "DelayedDataFrame", "DataFrame"; using the first of these

    gds_path <- tempfile(fileext=".gds")
    .initiate_snpgds(se, gds_path, compress="LZMA_RA")
    f <- gdsfmt::openfn.gds(gds_path)
    on.exit(gdsfmt::closefn.gds(f))
    expect_true(validObject(f))
    expect_true(is(f, "gds.class"))
    expect_true(all(c("sample.id", "snp.id", "snp.chromosome", "snp.position") %in% gdsfmt::ls.gdsn(f)))
})

test_that("write GDSSE works", {
    .initiate_seqgds <- VariantExperiment:::.initiate_seqgds
    .write_ve_as_gds <- VariantExperiment:::.write_ve_as_gds
    .write_ve_as_newve <- VariantExperiment:::.write_ve_as_newve
    file <- SeqArray::seqExampleFileName("gds")
    se <- makeVariantExperimentFromGDS(file, rowDataOnDisk=FALSE, colDataOnDisk=FALSE)
    gds_path <- tempfile(fileext=".gds")
    .initiate_seqgds(se, gds_path, compress="LZMA_RA")
    suppressWarnings(
        .write_ve_as_gds(se, "SEQ_ARRAY", gds_path, chunk_size = 10000,
                     compress = "LZMA_RA", verbose = FALSE))
    se1 <- .write_ve_as_newve(se, gds_path, "SEQ_ARRAY", TRUE, TRUE)
    expect_true(validObject(se1))
    expect_s4_class(se1, "VariantExperiment")
    expect_identical(dim(se), dim(se1))
})
