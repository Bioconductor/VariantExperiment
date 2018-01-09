file <- SeqArray::seqExampleFileName("gds")
se <- makeSummarizedExperimentFromGDS(file)
.get_gdsdata_fileFormat(gdsfile(se)[1])

gds_path <- tempfile(fileext=".gds")
gfile <- createfn.gds(filename=gds_path)
put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
put.attr.gdsn(gfile$root, "FileVersion", "v1.0")
addfolder.gdsn(gfile, "description", replace=TRUE)
add.gdsn(gfile, "sample.id", colnames(se))
add.gdsn(gfile, "variant.id", as.integer(rownames(se)), replace=T)
add.gdsn(gfile, "chromosome", as.character(seqnames(rowRanges(se))))
add.gdsn(gfile, "position", ranges(rowRanges(se))@start) ## ??
add.gdsn(gfile, "allele",
         paste(mcols(rowRanges(se))$REF,
               mcols(rowRanges(se))$ALT,
               sep=",")
         )  ## "," separated.

fill <- function(data, chunk_size = 10000) {
    path <- strsplit(data@seed@name, "/")[[1]]  ## gds node path.
    directories <- head(path, -1)
    rdnode <- index.gdsn(gfile, "")
    for (directory in directories) {
        if (!directory %in% ls.gdsn(index.gdsn(rdnode, ""))) {
            rdnode <- addfolder.gdsn(rdnode, name=directory, type="directory")
        } else {
            rdnode <- index.gdsn(rdnode, directory)
        }
    }
    ## FIXME: if a matrix / array, do this in 'chunks', e.g., of 10,000 rows
    idx <- seq(1, nrow(data), by = chunk_size)
    ridx = (idx[1] - 1) + seq_len( min(idx[1] + chunk_size, nrow(data) + 1) - idx[1] )
    if(length(dim(data)) == 1){
        value = unname(as.array(data[ridx, drop=FALSE]))
    }else if(length(dim(data) == 2)){
        value = unname(as.array(data[ridx,, drop=FALSE]))
    }else if(length(dim(data) == 3)){
        value = unname(as.array(data[ridx,,, drop=FALSE]))
    }
    ## need "drop"? if append later, no need for "drop=FALSE". 
    perm <- data@seed@permute
    if(perm) value <- aperm(value)
    itemnode <- add.gdsn(rdnode, tail(path, 1), val = value, compress = "LZMA_RA")
    if(length(idx) > 1){
        for (start in idx[-1]) {
            ridx = (start - 1) +
                seq_len( min(start + chunk_size, nrow(data) + 1) - start )
            if(length(dim(data)) == 1){
                value = unname(as.array(data[ridx, drop=FALSE]))
            }else if(length(dim(data)) == 2){
                value = unname(as.array(data[ridx,, drop=FALSE]))
            }else if(length(dim(data)) == 3){
                value = unname(as.array(data[ridx,,, drop=FALSE]))
            }
            if (perm) value = aperm(value)
            append.gdsn(itemnode, val=value, check=TRUE) ##?? FIXME?
        }
    }
}

for (elt in seq_along(assays(se)))
    fill(assays(se)[[elt]])
for (elt in names(rowData(se))[-c(2,3)])
    fill(rowData(se)[[elt]])
for (elt in names(colData(se)))
    fill(colData(se)[[elt]])

closefn.gds(gfile)
se1 <- makeSummarizedExperimentFromGDS(gds_path)
