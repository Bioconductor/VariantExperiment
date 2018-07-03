test_that("initiate seq gds file works", {
    .initiate_seqgds <- VariantExperiment:::.initiate_seqgds
    file <- SeqArray::seqExampleFileName("gds")
    se <- GDS2VE(file)
    gds_path <- tempfile(fileext=".gds")
    .initiate_seqgds(se, gds_path, compress="LZMA_RA")
    f <- gdsfmt::openfn.gds(gds_path)
    on.exit(gdsfmt::closefn.gds(f))
    expect_true(validObject(f))
    expect_true(class(f)== "gds.class")
    expect_true(all(c("sample.id", "variant.id", "chromosome", "position") %in% gdsfmt::ls.gdsn(f)))
})

test_that("initiate snp gds file works", {
    .initiate_snpgds <- VariantExperiment:::.initiate_snpgds
    file <- SNPRelate::snpgdsExampleFileName()
    se <- GDS2VE(file)
    ## FIXME: warning
    ##   1: In methods:::.selectDotsMethod(classes, .MTable, .AllMTable) :
    ## multiple direct matches: "DelayedDataFrame", "DataFrame"; using the first of these

    gds_path <- tempfile(fileext=".gds")
    .initiate_snpgds(se, gds_path, compress="LZMA_RA")
    f <- gdsfmt::openfn.gds(gds_path)
    on.exit(gdsfmt::closefn.gds(f))
    expect_true(validObject(f))
    expect_true(class(f)== "gds.class")
    expect_true(all(c("sample.id", "snp.id", "snp.chromosome", "snp.position") %in% gdsfmt::ls.gdsn(f)))
})

## test_that("write GDSSE works", {
##     .write_se_as_newse <- VariantExperiment::.write_se_as_newse
##     file <- SeqArray::seqExampleFileName("gds")
##     se <- GDS2VE(file, rowDataOnDisk=F, colDataOnDisk=F)
##     gds_path <- tempfile(fileext=".gds")
##     .initiate_seqgds(se, gds_path, compress="LZMA_RA")
##     se1 <- .write_se_as_newse(se, gds_path, "SEQ_ARRAY", TRUE, TRUE)
## })
