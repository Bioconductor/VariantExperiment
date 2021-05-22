.get_gds_fileFormat <- function(file)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    ff <- get.attr.gdsn(f$root)$FileFormat
    ff
}

.get_gdsnode_desp <- function(file, node, desp)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    objdesp <- objdesp.gdsn(index.gdsn(f, node))
    desp <- match.arg(desp, names(objdesp))
    objdesp[[desp]]
}

.get_gds_arraynodes <- function(gdsfile) {
    allnodes <- gdsnodes(gdsfile)
    isarray <- vapply(allnodes,
                      function(x) .get_gdsnode_desp(gdsfile, x, "is.array"),
                      logical(1))
    dim <- lapply(allnodes,
                  function(x) .get_gdsnode_desp(gdsfile, x, "dim"))
    res <- allnodes[isarray & lengths(dim) > 1 & 
                 ! vapply(dim, function(x) any(x == 0L), logical(1)) &
                 !grepl("~", allnodes)]
    res
}

## grep("chromosome", "position", "id") as
## GenomicRanges::Granges(seqnames=chromosome, ranges = Iranges(start
## = position, end = position + reflen-1L)) where define reflen = 1L
## for snps, but can be other values for other data type)

.get_gds_annonodes <- function(gdsfile, len.anno) {  
    allnodes <- gdsnodes(gdsfile)
    dim <- lapply(allnodes,
                  function(x) .get_gdsnode_desp(gdsfile, x, "dim"))
    idx <- lengths(dim) == 1 & ! vapply(dim, function(x) any(x == 0L), logical(1))
    res <- allnodes[idx][vapply(dim[idx], function(x) x[1] == len.anno, logical(1))]
    res ## returns character(0) if nothing matches
}

#' ShowAvailable
#' 
#' The function to show the available entries for the arguments within
#' \code{makeVariantExperimentFromGDS}
#' @name showAvailable
#' @rdname makeVariantExperimentFromGDS
#' @param file the path to the gds.class file.
#' @param feature.num the number of features for the gds file. Must be
#'     provided if the file format is not \code{SNP_ARRAY} or
#'     \code{SEQ_ARRAY}.
#' @param sample.num the number of features for the gds file. Must be
#'     provided if the file format is not \code{SNP_ARRAY} or
#'     \code{SEQ_ARRAY}.
#' @param args the arguments in \code{makeVariantExperimentFromGDS}.
#' @examples
#' ## snp gds file
#' gds <- SNPRelate::snpgdsExampleFileName()
#' showAvailable(gds)
#'
#' ## sequencing gds file
#' gds <- SeqArray::seqExampleFileName("gds")
#' showAvailable(gds)
#'
#' @importFrom IRanges CharacterList
#' @export
#' 
showAvailable <- function(file,  
                          args=c("assayNames", "rowDataColumns",
                                 "colDataColumns", "infoColumns"),
                          ftnode, smpnode)
{ 
    ## check if character.
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    args <- match.arg(args, several.ok=TRUE)
    ff <- .get_gds_fileFormat(file)
    ## if (is.null(ff))
    ##     stop("The gds file does not have the 'FileFormat' attribution!")
    if (!is.null(ff) && ff == "SEQ_ARRAY") {
        res <- .showAvailable_seqarray(file, args, ftnode = "variant.id", smpnode = "sample.id")
    } else if (!is.null(ff) && ff == "SNP_ARRAY") {
        res <- .showAvailable_snparray(file, args, ftnode = "snp.id", smpnode = "sample.id")
    } else {
        res <- .showAvailable_general(file, args, ftnode, smpnode)
    }
    res
}

.showAvailable_general <- function(file,
                                   args = c("assayNames",
                                            "rowDataColumns",
                                            "colDataColumns",
                                            "mcols"),
                                   ftnode, smpnode)
{
    res <- CharacterList()
    if("assayNames" %in% args) {
        res$assayNames <- .get_gds_arraynodes(file)
    }
    if("rowDataColumns" %in% args) {
        feature.num <- .get_gdsnode_desp(file, ftnode, "dim")
        res$rowDataColumns <- .get_gds_annonodes(file, len.anno = feature.num)
    }
    if ("colDataColumns" %in% args) {
        sample.num <- .get_gdsnode_desp(file, smpnode, "dim")
        res$colDataColumns <- .get_gds_annonodes(file, len.anno = sample.num)
    }
    res
}
