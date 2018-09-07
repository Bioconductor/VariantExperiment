## #' @rdname VariantExperiment-methods
## #' @export seqAlleleFreq seqAlleleCount seqMissing seqNumAllele

## SeqArray: seqAlleleFreq seqAlleleCount seqMissing seqNumAllele
setGeneric("seqAlleleFreq",
           function(gdsfile, ref.allele=0L, .progress=FALSE,
                    parallel=seqGetParallel())
               standardGeneric("seqAlleleFreq"),
           useAsDefault=SeqArray::seqAlleleFreq,
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
           useAsDefault=SeqArrayseqNumAllele,
           signature="gdsfile")
### note:
## seqNumAllele@generic
## [1] "seqNumAllele"
## attr(,"package")
## [1] "SeqArray"
## document()
## Creating a generic function for ‘seqNumAllele’ from package ‘SeqArray’ in package ‘VariantExperiment’

