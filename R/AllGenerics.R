## #' @export seqAlleleFreq seqAlleleCount seqMissing seqNumAllele hwe inbreedCoeff pca titv refDosage altDosage countSingletons homozygosity heterozygosity meanBySample missingGenotypeRate isSNV isVariant

## ## SeqArray: seqAlleleFreq seqAlleleCount seqMissing seqNumAllele
setGeneric("seqAlleleFreq",
           function(gdsfile, ref.allele=0L, .progress=FALSE,
                    parallel=seqGetParallel())
               standardGeneric("seqAlleleFreq"),
           useAsDefault=seqAlleleFreq,
           signature="gdsfile")
setGeneric("seqAlleleCount",
           function(gdsfile, ref.allele=0L, .progress=FALSE,
                    parallel=seqGetParallel())
               standardGeneric("seqAlleleCount"),
           useAsDefault=SeqArray::seqAlleleCount,
           signature="gdsfile")
setGeneric("seqMissing",
           function(gdsfile, per.variant=TRUE, .progress=FALSE,
                    parallel=seqGetParallel())
               standardGeneric("seqMissing"),
           useAsDefault=SeqArray::seqMissing,
           signature="gdsfile")
setGeneric("seqNumAllele",
           function(gdsfile)
               standardGeneric("seqNumAllele"),
           useAsDefault=seqNumAllele,
           signature="gdsfile")
### note:
## seqNumAllele@generic
## [1] "seqNumAllele"
## attr(,"package")
## [1] "SeqArray"
## document()
## Creating a generic function for ‘seqNumAllele’ from package ‘SeqArray’ in package ‘VariantExperiment’

## ## SeqVarTools: implemented: hwe, inbreedCoeff, pca, titv, refDosage, altDosage, countSingletons, homozygosity, heterozygosity, meanBySample, missingGenotypeRate, isSNV, isVariant
setGeneric("hwe",
           function(gdsobj, ...) standardGeneric("hwe"),
           useAsDefault=SeqVarTools::hwe,
           signature="gdsobj")
setGeneric("inbreedCoeff",
           function(gdsobj, ...) standardGeneric("inbreedCoeff"),
           useAsDefault=SeqVarTools::inbreedCoeff,
           signature="gdsobj")
setGeneric("pca",
           function(gdsobj, ...) standardGeneric("pca"),
           useAsDefault=SeqVarTools::pca,
           signature="gdsobj")
setGeneric("titv",
           function(gdsobj, ...) standardGeneric("titv"),
           useAsDefault=SeqVarTools::titv,
           signature="gdsobj")
setGeneric("refDosage",
           function(gdsobj, ...) standardGeneric("refDosage"),
           useAsDefault=SeqVarTools::refDosage,
           signature="gdsobj")
setGeneric("altDosage",
           function(gdsobj, ...) standardGeneric("altDosage"),
           useAsDefault=SeqVarTools::altDosage,
           signature="gdsobj")
setGeneric("countSingletons",
           function(gdsobj, ...) standardGeneric("countSingletons"),
           useAsDefault=SeqVarTools::countSingletons,
           signature="gdsobj")
setGeneric("homozygosity",
           function(gdsobj, ...) standardGeneric("homozygosity"),
           useAsDefault=SeqVarTools::homozygosity,
           signature="gdsobj")
setGeneric("heterozygosity",
           function(gdsobj, ...) standardGeneric("heterozygosity"),
           useAsDefault=SeqVarTools::heterozygosity,
           signature="gdsobj")
setGeneric("meanBySample",
           function(gdsobj, ...) standardGeneric("meanBySample"),
           useAsDefault=SeqVarTools::meanBySample,
           signature="gdsobj")
setGeneric("missingGenotypeRate",
           function(gdsobj, ...) standardGeneric("missingGenotypeRate"),
           useAsDefault=SeqVarTools::missingGenotypeRate,
           signature="gdsobj")
setGeneric("isSNV",
           function(gdsobj, ...) standardGeneric("isSNV"),
           useAsDefault=SeqVarTools::isSNV,
           signature="gdsobj")
setGeneric("isVariant",
           function(gdsobj, ...) standardGeneric("isVariant"),
           useAsDefault=SeqVarTools::isVariant,
           signature="gdsobj")

## 1 case: setGeneric:
## @importMethodsFrom seqArray
## @import SeqArray except=rowRanges
## .rowRanges <- function(x, ...){
##     seqArray::rowRanges(x)
## }
## setMethod("rowRanges", "VariantExperiment", .rowRanges)
## #' @export rowRanges
