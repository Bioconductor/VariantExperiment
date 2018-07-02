#' VariantExperiment: A package to represent VCF / GDS files using standard SummarizedExperiment metaphor with on-disk representation. 
#' @description The package \code{VariantExperiment} takes GDS file or VCF file as input, and save them in VariantExperiment object. Assay data are saved in \code{GDSArray} objects and annotation data are saved in \code{DelayedDataFrame} format, both of which remain on-disk until needed. Common manipulations like subsetting, mathematical transformation and statistical analysis are done easily and quickly in _R_. 
#' @docType package
#' @name VariantExperiment-package
#' @rdname VariantExperiment-package
NULL
