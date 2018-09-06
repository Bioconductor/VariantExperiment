# VariantExperiment

VariantExperiment is a Bioconductor package for saving data in VCF/GDS
format into RangedSummarizedExperiment object. The high-throughput
genetic/genomic data are saved in GDSArray objects. The annotation
data for features/samples are saved in DelayedDataFrame format with
mono-dimensional GDSArray in each column. The on-disk representation
of both assay data and annotation data achieves on-disk reading and
processing and saves memory space significantly. The interface of
RangedSummarizedExperiment data format enables easy and common
manipulations for high-throughput genetic/genomic data with common
SummarizedExperiment metaphor in R and Bioconductor.

