## #' GDSArraySeed
## #' Generate the seed for gds data for the GDSArray
## #' @import S4Vectors
## #' @import methods

### =========================================================================
### GDSArray objects
### -------------------------------------------------------------------------
## #' @importClassesFrom S4Vectors DataFrame
## ## #' @export 
## setClass("DelayedDataFrame",
##          contains = "DataFrame", ## from S4VectorsA data.frame-like interface for
##                                  ## S4 objects that implement length() and `[`
##          slots = c(
##              file="character",      ## Absolute path to the gds file so the object
##                                ## doesn't break when the user changes the working
##                                ## directory (e.g. with setwd()).
##              name="character", ## character vector, gds node names in the gds file.
##              dim = "integer",
##              dimnames = "list",  
##              permute = "logical",   ## logical vector.
##              first_val = "ANY"
##          )
##          )
## )



