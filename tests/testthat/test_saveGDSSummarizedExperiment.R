testthat("initiate seq gds file works", {
    .initiate_seqgds <- VariantExperiment::.initiate_seqgds
    file <- SeqArray::seqExampleFileName("gds")
    se <- makeSummarizedExperimentFromGDS(file)
    gds_path <- tempfile(fileext=".gds")
    .initiate_seqgds(se, gds_path, compress="LZMA_RA")
    f <- gdsfmt::openfn.gds(gds_path)
    expect_true(validObject(f))
    expect_true(class(f)== "gds.class")
    expect_true(all(c("sample.id", "variant.id", "chromosome", "position") %in% gdsfmt::ls.gdsn(f)))
    gdsfmt::closefn.gds(f)
})

testthat("initiate snp gds file works", {
    .initiate_snpgds <- VariantExperiment::.initiate_snpgds
    file <- SNPRelate::snpgdsExampleFileName()
    se <- makeSummarizedExperimentFromGDS(file)
    gds_path <- tempfile(fileext=".gds")
    .initiate_snpgds(se, gds_path, compress="LZMA_RA")
    f <- gdsfmt::openfn.gds(gds_path)
    expect_true(validObject(f))
    expect_true(class(f)== "gds.class")
    expect_true(all(c("sample.id", "snp.id", "snp.chromosome", "snp.position") %in% gdsfmt::ls.gdsn(f)))
    gdsfmt::closefn.gds(f)
})

## testthat("write GDSSE works", {
##     .write_se_as_newse <- VariantExperiment::.write_se_as_newse
##     file <- SeqArray::seqExampleFileName("gds")
##     se <- makeSummarizedExperimentFromGDS(file, rowDataOnDisk=F, colDataOnDisk=F)
##     gds_path <- tempfile(fileext=".gds")
##     .initiate_seqgds(se, gds_path, compress="LZMA_RA")
##     se1 <- .write_se_as_newse(se, gds_path, "SEQ_ARRAY", TRUE, TRUE)
## })
