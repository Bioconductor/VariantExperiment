## setMethod("rbind", "DelayedDataFrame", function (..., deparse.level = 1) 
## {
##     ## browser(); browser()
##     args <- list(...)
##     hasrows <- unlist(lapply(args, nrow), use.names = FALSE) > 
##         0L
##     hascols <- unlist(lapply(args, ncol), use.names = FALSE) > 
##         0L
##     if (!any(hasrows | hascols)) {
##         return(DataFrame())
##     }
##     else if (!any(hasrows)) {
##         return(args[[which(hascols)[1L]]])
##     }
##     else if (sum(hasrows) == 1) {
##         return(args[[which(hasrows)]])
##     }
##     else {
##         args <- args[hasrows]
##     }
##     df <- args[[1L]]
##     for (i in 2:length(args)) {
##         if (ncol(df) != ncol(args[[i]])) 
##             stop("number of columns for arg ", i, " do not match those of first arg")
##         if (!identical(colnames(df), colnames(args[[i]]))) 
##             stop("column names for arg ", i, " do not match those of first arg")
##     }
##     if (ncol(df) == 0) {
##         ans <- DataFrame()
##         ans@nrows <- sum(unlist(lapply(args, nrow), use.names = FALSE))
##     }
##     else {
##         cols <- lapply(colnames(df), function(cn) {
##             cols <- lapply(args, `[[`, cn)
##             isRle <- vapply(cols, is, logical(1L), "Rle")
##             if (any(isRle) && !all(isRle)) {
##                 cols[isRle] <- lapply(cols[isRle], decodeRle)
##             }
##             isFactor <- vapply(cols, is.factor, logical(1L))
##             if (any(isFactor)) {
##                 cols <- lapply(cols, as.factor)
##                 levs <- unique(unlist(lapply(cols, levels), use.names = FALSE))
##                 cols <- lapply(cols, factor, levs)
##             }
##             rectangular <- length(dim(cols[[1]])) == 2L
##             isDelayedArray <- vapply(cols, is, logical(1), "DelayedArray") ##
##             if (rectangular) {
##                 combined <- do.call(rbind, unname(cols))
##             } else {
##                 ## if DelayedArray, use "arbind" here?
##                 if (all(isDelayedArray)) {
##                     combined <- do.call("arbind", unname(cols))
##                 } else {
##                 combined <- do.call(c, unname(cols))
##                 }
##             }
##             if (any(isFactor)) 
##                 combined <- structure(combined, class = "factor", 
##                   levels = levs)
##             combined
##         })
##         names(cols) <- colnames(df)
##         ans <- S4Vectors:::new_DataFrame(cols, nrows = NROW(cols[[1]]))
##     }
##     rn <- unlist(lapply(args, rownames), use.names = FALSE)
##     if (!is.null(rn)) {
##         if (length(rn) != nrow(ans)) {
##             rn <- NULL
##         }
##         else if (anyDuplicated(rn)) 
##             rn <- make.unique(rn, sep = "")
##     }
##     rownames(ans) <- rn
##     if (!is.null(mcols(df))) {
##         df_mcols <- mcols(df)
##         if (all(sapply(args, function(x) identical(mcols(x), 
##             df_mcols)))) 
##             mcols(ans) <- df_mcols
##     }
##     ans
## })
