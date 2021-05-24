.permdim <- function(x, dim1, dim2) {
    pos <- match(c(dim1, dim2), dim(x))
    if (length(dim(x)) > 2) {
        aperm(x, perm = c(pos, setdiff(seq_along(dim(x)), pos)))
    } else {
        aperm(x, perm = pos)
    }
}

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

.get_gds_annonodes <- function(gdsfile, len.anno) {  
    allnodes <- gdsnodes(gdsfile)
    dim <- lapply(allnodes,
                  function(x) .get_gdsnode_desp(gdsfile, x, "dim"))
    idx <- lengths(dim) == 1 & ! vapply(dim, function(x) any(x == 0L), logical(1))
    res <- allnodes[idx][vapply(dim[idx], function(x) x[1] == len.anno, logical(1))]
    res ## returns character(0) if nothing matches
}
