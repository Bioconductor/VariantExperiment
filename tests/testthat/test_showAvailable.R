snpfile <- SNPRelate::snpgdsExampleFileName()
seqfile <- SeqArray::seqExampleFileName("gds")
gdsfile <- system.file("extdata/test.gds", package = "VariantExperiment")

test_that("showAvailable for snp gds file works", {
    .showAvailable_snparray <- VariantExperiment:::.showAvailable_snparray
    res <- .showAvailable_snparray(snpfile)
    expect_equal(length(res), 3L)
    expect_equal(res$assayNames, "genotype")
    expect_equal(res$rowDataColumns, c("snp.rs.id", "snp.allele"))
    expect_equal(res$colDataColumns, c("family.id", "father.id", "mother.id", "sex", "pop.group"))
})

test_that("showAvailable for seq gds file works", {
    .showAvailable_seqarray <- VariantExperiment:::.showAvailable_seqarray
    res <- .showAvailable_seqarray(seqfile)
    expect_equal(length(res), 4L)
    expect_equal(names(res), c("assayNames", "rowDataColumns", "colDataColumns", "infoColumns"))
    expect_equal(res$assayNames, c("genotype/data", "phase/data", "annotation/format/DP/data"))
    expect_equal(res$rowDataColumns, c("allele", "annotation/id", "annotation/qual", "annotation/filter"))
    expect_equal(res$colDataColumns, "family")
    expect_equal(res$infoColumns, c("AC", "AN", "DP", "HM2", "HM3", "OR", "GP", "BN"))
})

test_that("showAvailable for general gds file works", {
    .showAvailable_general <- VariantExperiment:::.showAvailable_general
    res <- .showAvailable_general(gdsfile, ftnode = "feature.id", smpnode = "sample.id")
    expect_equal(length(res), 3L)
    expect_equal(res$assayNames, "matrix")
    expect_equal(res$rowDataColumns, c("feature.id", "chromosome", "position"))
    expect_equal(res$colDataColumns, c("sample.id", "sample.annotation/family"))
})

