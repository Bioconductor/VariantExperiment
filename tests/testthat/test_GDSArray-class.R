test_that("GDSArray constructor works", {
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    gds <- GDSArray(file)
    expect_s4_class(gds, "GDSArray")
    expect_true(validObject(gds))
    expect_equal(dim(gds), c(9088L, 279L))
})

test_that("info input works", {
    .info_seqarray_ondisk <- VariantExperiment:::.info_seqarray_ondisk
    file <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")

    df <- .info_seqarray_ondisk(file, character())
    expect_equal(dim(df), c(1348L, 9L))
    expect_equal(names(df)[[1]], "info_AA")
    expect_equal(class(df[["info_AA"]][[1]]), "character")
})
