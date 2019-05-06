test_that("Allele related functions work", {
    vcf <- SeqArray::seqExampleFileName("vcf")
    sample.info <- system.file("extdata", "Example_sampleInfo.txt",
                               package="VariantExperiment")
    ve <- suppressWarnings(
        makeSummarizedExperimentFromVCF(vcf, out.dir = tempfile(),
                                        sample.info = sample.info))
        
        ## seqAlleleFreq
        freqAll <- seqAlleleFreq(ve, ref.allele=NULL)
        expect_true(validObject(freqAll))
        expect_equal(length(freqAll), nrow(ve))
        expect_equal(lengths(freqAll)[1338], 3)
        expect_equal(lengths(freqAll)[1323], 3)

        freq.ref <- seqAlleleFreq(ve, ref.allele = 0L)
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
        expect_equal(nm.allele,
                     c(rep(2L, 1322), 3L, rep(2L, 14), 3L, rep(2L, 10)))

        ## seqMissing
        mr.var <- seqMissing(ve, per.variant=TRUE)
        expect_equal(length(mr.var), nrow(ve))

        mr.samp <- seqMissing(ve, per.variant=FALSE)
        expect_equal(length(mr.samp), ncol(ve))

    })

### hwe, inbreedCoeff, pca, titv, refDosage, altDosage,
### countSingletons, heterozygosity, homozygosity, meanBySample,
### missingGenotypeRate (by sample/variant), isSNV, isVariant

test_that("other statistical functions work", {
    vcf <- SeqArray::seqExampleFileName("vcf")
    sample.info <- system.file("extdata", "Example_sampleInfo.txt",
                               package="VariantExperiment")
    ve <- suppressWarnings(
        makeSummarizedExperimentFromVCF(vcf, out.dir = tempfile(),
                                        sample.info = sample.info))

        ## hwe
        h <- hwe(ve)
        expect_equal(dim(h), c(1348L, 7L))
        expect_equal(colnames(h),
                     c("variant.id", "nAA", "nAa", "naa", "afreq", "p", "f"))

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
        expect_identical(dim(dosage.ref), dim(ve))
        expect_equal(rownames(dosage.ref), rownames(ve))
        expect_equal(colnames(dosage.ref), colnames(ve))

        ## altDosage
        dosage.alt <- altDosage(ve)
        expect_identical(dim(dosage.alt), dim(ve))
        dosage.alt <- altDosage(ve, sparse=TRUE)
        expect_identical(dim(dosage.alt), dim(ve))    
        expect_s4_class(dosage.alt, "dgCMatrix")
        
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
        mn <- meanBySample(ve, var.name="annotation/format/DP/data",
                           use.names=TRUE)
        expect_equal(names(mn), colnames(ve))
        expect_error(meanBySample(ve, var.name="genotype", use.names=TRUE))
        
        ## isSNV
        issnv <- isSNV(ve, biallelic=TRUE)    ## Setting ‘biallelic=TRUE’
        ## is considerably faster
        ## for large data sets.
        expect_true(is.logical(issnv))
        expect_equal(length(issnv), nrow(ve))
        
        ## isVariant
        isvar <- isVariant(ve)
        expect_true(is.logical(isvar))
        expect_identical(dim(isvar), dim(ve))
    })

