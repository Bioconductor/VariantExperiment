test_that("VCF2VE works", {
    vcf <- SeqArray::seqExampleFileName("vcf")
    ve <- VCF2VE(vcf, out.dir = tempfile())

    ## general
    expect_equal(dim(ve), c(1348L, 90L))
    expect_equal(dim(rowData(ve)), c(1348L, 14L))
    expect_equal(dim(colData(ve)), c(90L, 0L))
    expect_s4_class(rowData(ve), "DelayedDataFrame")
    expect_s4_class(colData(ve), "DelayedDataFrame")
    expect_equal(rownames(rowData(ve)), rownames(ve))
    expect_equal(rownames(colData(ve)), colnames(ve))
    
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

    ve2 <- VCF2VE(vcf, out.dir = tempfile(), fmt.import="xx")
    expect_equivalent(ve, ve2)
    ve3 <- VCF2VE(vcf, out.dir = tempfile(), fmt.import=NULL)
    expect_equivalent(ve, ve3)
    ve4 <- VCF2VE(vcf, out.dir = tempfile(), fmt.import="DP")
    expect_equivalent(ve, ve4)
    ve6 <- VCF2VE(vcf, out.dir = tempfile(), reference="hg19")
    expect_equivalent(ve, ve6)
    
    ve7 <- VCF2VE(vcf, out.dir = tempfile(), start=101, count=1000) 
    expect_equal(dim(ve7), c(1000L, 90L))
    expect_equal(rownames(ve7), as.character(seq(101, length.out=1000)))


    sample.info <- system.file("extdata", "Example_sampleInfo.txt", package="VariantExperiment")
    ve8 <- VCF2VE(vcf, out.dir = tempfile(), sample.info = sample.info)
    expect_equal(dim(colData(ve8)), c(90L, 1L))
    expect_equal(colnames(colData(ve8)), "family")
})

test_that("Allele related functions work", {
    vcf <- SeqArray::seqExampleFileName("vcf")
    sample.info <- system.file("extdata", "Example_sampleInfo.txt", package="VariantExperiment")
    ve <- VCF2VE(vcf, out.dir = tempfile(), sample.info = sample.info)

    ## seqAlleleFreq
    freqAll <- seqAlleleFreq(ve, ref.allele=NULL)
    expect_true(validObject(freqAll))
    expect_equal(length(freqAll), nrow(ve))
    expect_equal(lengths(freqAll)[1338], 3)
    expect_equal(lengths(freqAll)[1323], 3)

    freq.ref <- seqAlleleFreq(ve)
    ref.allele <- as.array(rowData(ve)$REF)
    freq.ref1 <- seqAlleleFreq(ve, ref.allele = ref.allele)
    expect_equal(freq.ref, freq.ref1)

    freq.alt <- seqAlleleFreq(ve, ref.allele = 1L)
    alt.allele <- sub(",[ATCG]*", "", as.array(rowData(ve)$ALT))
    freq.alt1 <- seqAlleleFreq(ve, ref.allele = alt.allele)
    expect_equal(freq.alt, freq.alt1)

    ## seqAlleleCount
    ct.ref <- seqAlleleCount(ve, ref.allele = 0L)
    ct.ref1 <- seqAlleleCount(ve, ref.allele = ref.allele)
    expect_equal(ct.ref, ct.ref1)

    ct.alt <- seqAlleleCount(ve, ref.allele = 1L)
    ct.alt1 <- seqAlleleCount(ve, ref.allele = alt.allele)
    expect_equal(ct.alt, ct.alt1)

    ## seqNumAllele
    nm.allele <- seqNumAllele(ve)
    expect_equal(nm.allele, c(rep(2L, 1322), 3L, rep(2L, 14), 3L, rep(2L, 10)))

    ## seqMissing
    mr.var <- seqMissing(ve, per.variant=TRUE)
    expect_equal(length(mr.var), nrow(ve))

    mr.samp <- seqMissing(ve, per.variant=FALSE)
    expect_equal(length(mr.samp), ncol(ve))

})

### hwe, inbreedCoeff, pca, titv, refDosage, altDosage, countSingletons, heterozygosity, homozygosity, meanBySample, missingGenotypeRate (by sample/variant), isSNV, isVariant

test_that("hwe works", {
    vcf <- SeqArray::seqExampleFileName("vcf")
    sample.info <- system.file("extdata", "Example_sampleInfo.txt", package="VariantExperiment")
    ve <- VCF2VE(vcf, out.dir = tempfile(), sample.info = sample.info)

    ## hwe
    h <- hwe(ve)
    expect_equal(dim(h), c(1348L, 7L))
    expect_equal(colnames(h), c("variant.id", "nAA", "nAa", "naa", "afreq", "p", "f"))

    ## inbreedCoeff
    inb.samp <- inbreedCoeff(ve, margin="by.sample", use.names=TRUE)
    expect_equal(names(inb.samp), colnames(ve))
    inb.var <- inbreedCoeff(ve, margin="by.variant")
    expect_equal(length(inb.var), nrow(ve))

    ## pca
    p <- pca(ve)
    expect_equal(dim(p$eigenvect), c(90L, 32L))
    expect_equal(length(p$eigenval), 32L)
    expect_equal(rownames(p$eigenvect), colnames(ve))

    ## titv
    titv.samp <- titv(ve, by.sample=TRUE, use.names=TRUE)
    expect_equal(names(titv.samp), colnames(ve))
    titv.all <- titv(ve, by.sample = FALSE)
    expect_equal(length(titv.all), 1L)

    ## refDosage
    dosage.ref <- refDosage(ve)
    expect_equal(dim(dosage.ref), c(90L, 1348L))
    expect_equal(rownames(dosage.ref), colnames(ve))
    expect_equal(colnames(dosage.ref), rownames(ve))

    ## altDosage
    dosage.alt <- altDosage(ve)  ## FIXME: sparse=TRUE does not work...

    ## coutnSingletons
    cts <- countSingletons(ve) 
    expect_equal(length(cts), ncol(ve))

    ## heterozygosity
    hetr <- heterozygosity(ve, margin="by.variant")
    expect_equal(length(hetr), nrow(ve))

    hetr <- heterozygosity(ve, margin="by.sample")
    expect_equal(length(hetr), ncol(ve))

    ## homozygosity
    homr.ref <- homozygosity(ve, allele="ref", margin="by.variant")
    expect_equal(length(homr.ref), nrow(ve))

    homr.alt <- homozygosity(ve, allele="alt", margin="by.variant")
    homr.all <- homozygosity(ve, allele="any", margin="by.variant")
    expect_true(max(homr.alt) == 1)
    expect_true(max(homr.all) == 1)
    
    ## meanBySample
    mn <- meanBySample(ve, var.name="annotation/format/DP", use.names=TRUE)
    ## FIXME: "var.name" need to be consistent with VE infrastructure. 
    expect_equal(names(mn), colnames(ve))

    ## missingGenotypeRate
    mr <- missingGenotypeRate(ve, margin="by.sample")
    expect_identical(mr, seqMissing(ve, per.variant=FALSE))
    
    mr <- missingGenotypeRate(ve, margin="by.variant")
    expect_identical(mr, seqMissing(ve, per.variant=TRUE))

    ## isSNV
    issnv <- isSNV(ve, biallelic=TRUE)    ## Setting ‘biallelic=TRUE’
                                          ## is considerably faster
                                          ## for large data sets.
    expect_true(is.logical(issnv))
    expect_equal(length(issnv), nrow(ve))
    
    ## isVariant
    isvar <- isVariant(ve)
    expect_true(is.logical(isvar))
    expect_equal(dim(isvar), c(90L, 1348L))
})
