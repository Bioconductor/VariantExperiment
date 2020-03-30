test_that("info input works", {
    .info_seqgds <- VariantExperiment:::.info_seqgds
    file##  <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    ## df <- .info_seqgds(file, character(), TRUE)
    ## expect_equal(dim(df), c(1348L, 9L))
    ## expect_equal(names(df)[[1]], "info_AA")
    ## expect_equal(class(df[["info_AA"]][[1]]), "character")

    ## df <- .info_seqgds(file, character(), FALSE)
    ## expect_s4_class(df[["info_AA"]], "CharacterList")

    ## expect_warning(.info_seqgds(file, "random", TRUE))
    ## expect_warning(.info_seqgds(file, "random", FALSE))
})

test_that("grange generation works", {
    .granges_gdsdata <- VariantExperiment:::.granges_gdsdata
    file <- SNPRelate::snpgdsExampleFileName()
    gr <- .granges_gdsdata(file, "SNP_ARRAY")
    expect_s4_class(gr, "GRanges")
    expect_s4_class(ranges(gr), "IRanges")
})

test_that("alleles for snpgds works", {
    .varnode_snpgds_inmem <- VariantExperiment:::.varnode_snpgds_inmem
    file <- SNPRelate::snpgdsExampleFileName()
    als <- .varnode_snpgds_inmem(file, "allele")
    expect_s4_class(als$ALLELE1, "DNAStringSet")
    expect_s4_class(als$ALLELE2, "DNAStringSet")
    expect_equal(dim(als), c(9088L, 2L))
})

test_that("alleles for seqgds works", {
    .varnode_seqgds_inmem <- VariantExperiment:::.varnode_seqgds_inmem
    file <- SeqArray::seqExampleFileName("gds")
    expect_equal(class(.varnode_seqgds_inmem(file, "id")[[1]]), "character")
    expect_s4_class(.varnode_seqgds_inmem(file, "ref")[[1]], "DNAStringSet")
    expect_s4_class(.varnode_seqgds_inmem(file, "alt")[[1]], "DNAStringSetList")
    expect_equal(class(.varnode_seqgds_inmem(file, "qual")[[1]]), "numeric")
    expect_equal(class(.varnode_seqgds_inmem(file, "filter")[[1]]), "character")
})

test_that("variant related nodes on disk reading works", {
    .varnode_gdsdata_ondisk <- VariantExperiment:::.varnode_gdsdata_ondisk
    file <- SeqArray::seqExampleFileName("gds")
    id <- .varnode_gdsdata_ondisk(file, "SEQ_ARRAY", "id")
    expect_s4_class(id, "GDSArray")
})

test_that("sample related nodes on disk reading works", {
    .sampnode_gdsdata_ondisk <- VariantExperiment:::.sampnode_gdsdata_ondisk
    file <- SeqArray::seqExampleFileName("gds")
    annot <- .sampnode_gdsdata_ondisk(file, "SEQ_ARRAY", "family")
    expect_s4_class(annot, "GDSArray")
})

test_that("showing available arguments works", {
   file <- SeqArray::seqExampleFileName("gds")
   expect_true(is(showAvailable(file), "CharacterList"))
   file1 <- SNPRelate::snpgdsExampleFileName()
   expect_null(showAvailable(file1)$infoColumns)
})
