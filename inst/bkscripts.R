.info_seqarray <- function(seqArrayFile, info){
    f <- SeqArray::seqOpen(seqArrayFile)
    on.exit(SeqArray::seqClose(f))
    infodata <- SeqArray::info(f)
    infonodes <- gdsfmt::ls.gdsn(index.gdsn(f, "annotation/info"))
    if(length(info)>0){
        idx <- info %in% infonodes
        if(any(!idx)){
            warning("\n", 'The "info" argument of "',
                    paste(info[!idx], collapse = ", "),
                    '" does not exist!', "\n",
                    'Please use showAvailable(file, "info") ',
                    'to get the available columns for "info."', "\n")
        }
        infodata <- infodata[names(infodata) %in% info]
    }
    infodata
}

.fmt_seqarray <- function(seqArrayFile, fmt, fmtOnDisk){
    f <- SeqArray::seqOpen(seqArrayFile)
    on.exit(SeqArray::seqClose(f))
    fmtnodes <- gdsfmt::ls.gdsn(index.gdsn(f, "annotation/format"))
    res <- seqGetData(f, "annotation/format/DP")  ## any other format nodes except DP?
    res$length <- as.matrix(res$length)
    res$data <- t(res$data)
    if(fmtOnDisk){
        res$length <- DelayedArray(res$length)
        res$data <- DelayedArray(res$data)
    }
    if(length(fmt)>0){
        idx <- fmt %in% fmtnodes
        if(any(!idx)){
            warning("\n", 'The "fmt" argument of "',
                    paste(fmt[!idx], collapse = ", "),
                    '" does not exist!', "\n",
                    'Please use showAvailable(file, "fmt") ',
                    'to get the available columns for "fmt."', "\n")
        }
    }
    res      
}


#' @importFrom Biostrings DNAStringSet
#' @import SNPRelate
#' @import SeqArray
.rowRanges_gdsdata <- function(file, fileFormat, rowDataColumns){
    if(fileFormat == "SNP_ARRAY"){
        f <- SNPRelate::snpgdsOpen(file)
        on.exit(SNPRelate::snpgdsClose(f))
    }else if(fileFormat == "SEQ_ARRAY"){
        f <- SeqArray::seqOpen(file)
        on.exit(SeqArray::seqClose(f))
    }
    rowDataColumns <- toupper(rowDataColumns)
    ## if(length(rowDataColumns) == 0){
    ##     idx.var <- setNames(rep(TRUE, 6), c("id", "allele", "ref", "alt", "qual", "filter")) 
    ## }else{
    ##     idx.var <- c(
    ##         id = any(grepl("ID", rowDataColumns)),
    ##         allele = any(grepl("ALLELE", rowDataColumns)),
    ##         ref = any(grepl("REF", rowDataColumns)),
    ##         alt = any(grepl("ALT", rowDataColumns)),
    ##         qual = any(grepl("QUAL", rowDataColumns)),
    ##         filter = any(grepl("FILTER", rowDataColumns))
    ##     )
    ## }  ## ? by default, read in all. 
    if(inherits(f, "SNPGDSFileClass")){
        meta <- DataFrame(
            ID = read.gdsn(index.gdsn(f, "snp.rs.id")),
            ALLELE1 = Biostrings::DNAStringSet(.alleles_snpgds(f)$allele1),
            ALLELE2 = Biostrings::DNAStringSet(.alleles_snpgds(f)$allele2)
            ## DNAStringSetList class
        )
        ## idx.var <- idx.var[1:2]
        ## if(any(!idx.var)){
        ##     rm.var <- names(idx.var[!idx.var])
        ##     meta <- meta[-match(toupper(rm.var), names(meta))]
        ## }
        idx.other <- rowDataColumns %in% c("ID", "ALLELE")
        if(any(!idx.other)){
            warning('The snp annotation of "',
                    tolower(paste(rowDataColumns[!idx.other], collapse = ", ")),
                    '" does not exist!',
                    'Please use showAvailable(file, "colDataColumns") ',
                    'to get the available columns for "colData."', "\n")
        }
        rr <- .granges_snpgds(f)
        mcols(rr) <- meta
    }else if(inherits(f, "SeqVarGDSClass")){
        meta <- DataFrame(ID = SeqArray::seqGetData(f, "annotation/id"),
                          REF = SeqArray::ref(f),
                          ALT = SeqArray::alt(f),
                          QUAL = SeqArray::qual(f),
                          FILTER = SeqArray::filt(f),
                          SeqArray::info(f))
        
        ##
        ## info <- SeqArray::info(f), only read querable in "makeSEFromGDS(info=)".
        ## Save all others as list elements. DataFrame(list).
        ## GDSArray(GDSArraySeed) for ref(f), alt(f), qual(f), filt(f).
        ## direct conversion of info(f) as DataFrame. 
        
        ## idx.var <- idx.var[-2]
        ## if(any(!idx.var)){
        ##     rm.var <- names(idx.var[!idx.var])
        ##     meta <- meta[-match(toupper(rm.var), names(meta))]
        ## }
        idx.other <- rowDataColumns %in% c("ID", "REF", "ALT", "QUAL", "FILTER")
        if(any(!idx.other)){
            warning('The variant annotation of "',
                    tolower(paste(rowDataColumns[!idx.other], collapse = ", ")),
                    '" does not exist!',
                    'Please use showAvailable(file, "rowDataColumns") ',
                    'to get the available columns for "rowData."', "\n")
        }
        rr <- SeqArray::granges(f)
        mcols(rr) <- meta
    }
    rr
}    
