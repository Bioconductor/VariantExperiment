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

#' ShowAvailable
#' 
#' The function to show the available entries for the arguments within
#' \code{makeVariantExperimentFromGDS}
#' @name showAvailable
#' @param file the path to the gds.class file.
#' @param args the arguments in \code{makeVariantExperimentFromGDS}.
#' @param ftnode the node name for feature id (e.g., "variant.id",
#'     "snp.id", etc.). Must be provided if the file format is not
#'     \code{SNP_ARRAY} or \code{SEQ_ARRAY}.
#' @param smpnode the node name for sample id (e.g.,
#'     "sample.id"). Must be provided if the file format is not
#'     \code{SNP_ARRAY} or \code{SEQ_ARRAY}.
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
