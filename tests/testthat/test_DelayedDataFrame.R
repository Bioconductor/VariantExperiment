context("DelayedDataFrame")

test_that("DelayedDataFrame constructor works", {
    obj <- DelayedDataFrame()
    expect_true(validObject(obj))
    expect_identical(c(0L, 0L), dim(obj))

    obj <- DelayedDataFrame(letters, LETTERS)
    expect_true(validObject(obj))
    expect_identical(c(26L, 2L), dim(obj))
    expect_identical(list(NULL, c("letters", "LETTERS")), dimnames(obj))

    obj <- DelayedDataFrame(letters, LETTERS, row.names=LETTERS)
    expect_identical(LETTERS, rownames(obj))

    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(letters, da1=I(da1), da2=I(da2), da3 = I(da3))
    expect_true(validObject(obj))
    expect_identical(
        list(NULL, c("letters", "da1", "da2", "da3")),
        dimnames(obj)
    )
    exp <- list(da1 = c(26L, 1L), da2 = c(26L, 2L), da3 = c(26L, 1L, 1L))
    expect_identical(exp, lapply(obj[-1], dim))
})

test_that("DelayedDataFrame [[ works", {
    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(letters, da1=I(da1), da2=I(da2), da3 = I(da3))

    expect_identical(letters, obj[["letters"]])
    expect_equal(da1, obj[["da1"]])     # FAILS: 'package' attribute stripped
    expect_equivalent(da2, obj[["da2"]])
    expect_equivalent(da3, obj[["da3"]])

    expect_identical(letters, obj[[1]])
    expect_equivalent(da1, obj[[2]])
    expect_equivalent(da2, obj[[3]])
    expect_equivalent(da3, obj[[4]])

    expect_error(obj[[5]], "subscript is out of bounds")
    expect_identical(NULL, obj[["da4"]]) # data.frame() behaves this way
    expect_error(obj[[-1]], "\\[\\[ subscript must be >= 1")
})

test_that("DelayedDataFrame 1-dimensional [ works", {
    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(letters, da1=I(da1), da2=I(da2), da3 = I(da3))

    ## empty
    expect_identical(obj, obj[])

    ## numeric
    expect_identical(DelayedDataFrame(letters), obj[1])
    expect_identical(DelayedDataFrame(da1 = I(da1)), obj[2])
    expect_identical(DelayedDataFrame(letters, da1 = I(da1)), obj[1:2])
    expect_identical(DelayedDataFrame(da1 = I(da1), letters), obj[2:1])
    expect_identical(DelayedDataFrame(letters, da1 = I(da1)), obj[-(3:4)])

    ## character
    expect_identical(DelayedDataFrame(letters), obj["letters"])
    expect_identical(DelayedDataFrame(da1 = I(da1)), obj["da1"])
    expect_identical(
        DelayedDataFrame(letters, da1 = I(da1)),
        obj[c("letters", "da1")]
    )
    expect_identical(
        DelayedDataFrame(da1 = I(da1), letters),
        obj[c("da1", "letters")]
    )

    ## logical
    expect_identical(
        DelayedDataFrame(letters, da2 = I(da2)),
        obj[c(TRUE, FALSE)]
    )
})

test_that("DelayedDataFrame 2-dimensional [ works", {
    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(letters, da1=I(da1), da2=I(da2), da3 = I(da3))

    ## empty
    expect_identical(obj, obj[,])

    ## columns
    expect_identical(letters, obj[,1])
    expect_identical(DelayedDataFrame(letters), obj[, 1, drop=FALSE])
    
    expect_equivalent(da1, obj[,2])
    expect_equivalent(da2, obj[,3])
    expect_equivalent(da3, obj[,4])

    expect_equivalent(DelayedDataFrame(da1 = I(da1)), obj[,2, drop=FALSE])
    expect_equivalent(DelayedDataFrame(da2 = I(da2)), obj[,3, drop=FALSE])
    expect_equivalent(DelayedDataFrame(da3 = I(da3)), obj[,4, drop=FALSE])

    expect_equivalent(DelayedDataFrame(letters, da1 = I(da1)), obj[,1:2])
    expect_equivalent(DelayedDataFrame(da1 = I(da1), letters), obj[,2:1])

    ## column names
    expect_equivalent(
        DelayedDataFrame(letters, da1 = I(da1)),
        obj[,c("letters", "da1")]
    )

    ## logical
    expect_identical(
        DelayedDataFrame(letters, da2 = I(da2)),
        obj[,c(TRUE, FALSE)]
    )
    
    ## rows
    ddf <- function(i)
        DelayedDataFrame(
            letters = letters[i],
            da1 = I(da1[i,]), da2 = I(da2[i,]), da3 = I(da3[i,,])
        )
    idx <- 1
    expect_equivalent(ddf(idx), obj[idx,])
    idx <- 1:5
    expect_equivalent(ddf(idx), obj[idx,])
    idx <- sample(nrow(obj))
    expect_equivalent(ddf(idx), obj[idx,])

    ## rownames
    obj <- DelayedDataFrame(
        letters, da1=I(da1), da2=I(da2), da3 = I(da3),
        row.names = letters
    )
    idx <- sample(rownames(obj))
    expect_equivalent(ddf(match(idx, rownames(obj))), obj[idx,])

    ## logical
    idx <- c(TRUE, FALSE)
    expect_equivalent(ddf(idx), obj[idx,])
})

test_that("[<-,DelayedDataFrame 2-D subseting works", {
    da0 <- DelayedArray(array(1:26, 26))
    obj <- DelayedDataFrame(letters, da0 = I(da0))

    obj[1:5, 1] <- rev(letters[1:5])
    expect_identical(rev(letters[1:5]), obj[[1]][1:5])

    ## TODO: subset-replace a DelayedArray column
})

test_that("cbind,DelayedDataFrame works", {
    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(letters, da1=I(da1), da2=I(da2), da3 = I(da3))

    expect_equivalent(obj, cbind(obj[1:2], obj[3:4]))
})

test_that("rbind,DelayedDataFrame works", {
    da0 <- DelayedArray(array(1:26, 26))
    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(
        letters, da0 = I(da0), da1=I(da1), da2=I(da2), da3 = I(da3)
    )

    expect_equivalent(obj, rbind(obj[1:10,], obj[11:26,]))
})

test_that("arbind,DelayedDataFrame works", {
    da0 <- DelayedArray(array(1:26, 26))
    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(
        letters, da0 = I(da0), da1=I(da1), da2=I(da2), da3 = I(da3)
    )

    expect_equivalent(obj, arbind(obj[1:10,], obj[11:26,]))
})

test_that("acbind,DelayedDataFrame works", {
    da0 <- DelayedArray(array(1:26, 26))
    da1 <- DelayedArray(matrix(1:26, 26, 1))
    da2 <- DelayedArray(matrix(1:52, 26, 2))
    da3 <- DelayedArray(array(1:26, c(26, 1, 1)))
    obj <- DelayedDataFrame(
        letters, da0 = I(da0), da1=I(da1), da2=I(da2), da3 = I(da3)
    )

    expect_equivalent(obj, acbind(obj[1:10,], obj[11:26,]))
})

context("LazyIndex")

test_that(".lazyIndex_inuse works", {
    ## 0. No updates
    ll <- .LazyIndex(list(1:10, 10:1, 11:20), index=1:3)
    expect_identical(ll, .lazyIndex_inuse(ll))

    ## 1. duplicate listData
    ll <- .LazyIndex(list(1:10, 10:1, 1:10), index=1:3)
    exp <- .LazyIndex(list(1:10, 10:1), index = c(1L, 2L, 1L))
    expect_identical(exp, .lazyIndex_inuse(ll))
    
    ll <- .LazyIndex(list(1:10, 10:1, 1:10), index=c(1:3, 3L))
    exp <- .LazyIndex(list(1:10, 10:1), index = c(1L, 2L, 1L, 1L))
    expect_identical(exp, .lazyIndex_inuse(ll))

    ## 2. listData in not in use index
    ll <- .LazyIndex(list(1:10, 10:1, 11:20), index=c(1L, 3L, 3L))
    exp <- .LazyIndex(list(1:10, 11:20), index = c(1L, 2L, 2L))
    expect_identical(exp, .lazyIndex_inuse(ll))
    
    ## 3. reorder
    ll <- .LazyIndex(list(1:10, 10:1, 11:20), index=c(3L, 2L, 1L))
    exp <- .LazyIndex(list(11:20, 10:1, 1:10), index=1:3)
    expect_identical(exp, .lazyIndex_inuse(ll))
})

test_that(".update_index works", {
    ll <- .LazyIndex(list(1:10, 10:1), index=1:2)

    lazyIndex <- sample(1:10)
    exp <- .LazyIndex(list(lazyIndex, 10:1), index=1:2)
    expect_identical(exp, .update_index(ll, 1, lazyIndex))

    exp <- .LazyIndex(list(10:1), index=c(1L, 1L))
    expect_identical(exp, .update_index(ll, 1, 10:1))

    exp <- .LazyIndex(list(NULL, 10:1), index=1:2)
    expect_identical(exp, .update_index(ll, 1, NULL))
})

test_that("LazyIndex [<- 1-D accounting is correct", {
    da0 <- DelayedArray(array(1:26, 26))
    obj <- DelayedDataFrame(letters, da0 = I(da0))

    obj <- obj[10:1,]
    obj["da0"] <- da0[sample(10)]
    exp <- .LazyIndex(list(10:1, NULL), index=1:2)
    expect_identical(exp, obj@lazyIndex)
})

test_that("LazyIndex [<- 2-D accounting is correct", {
    da0 <- DelayedArray(array(1:26, 26))
    obj <- DelayedDataFrame(letters, da0 = I(da0))

    obj[1:5, 1] <- rev(letters[1:5])
    nullList <- .LazyIndex(vector("list", 1), index=rep(1L, length(obj)))
    expect_identical(nullList, obj@lazyIndex)

    da0 <- DelayedArray(array(1:26, 26))
    obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj <- obj[10:1,]
    obj[1:5, 1] <- rev(letters[1:5])
    exp <- .LazyIndex(list(NULL, 10:1), index=1:2)
    expect_identical(exp, obj@lazyIndex)
    expect_identical(letters[c(5:1, 5:1)], obj@listData[[1]])
})

test_that("lazyIndex [[<- accounting is correct", {
    da0 <- DelayedArray(array(1:26, 26))
    obj <- DelayedDataFrame(letters, da0 = I(da0))

    nullList <- .LazyIndex(vector("list", 1), index=rep(1L, length(obj)))
    expect_identical(nullList, obj@lazyIndex)

    obj <- obj[1:10,]
    exp <- .LazyIndex(list(1:10), index=rep(1L, length(obj)))
    expect_identical(exp, obj@lazyIndex)

    obj0 <- obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj <- obj[1:10,]
    expect_identical(nullList, obj0@lazyIndex)

    ## .append_list_element
    obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj0 <- obj <- obj[1:10,]
    obj[["x"]] <- da0[1:10]
    exp <- .LazyIndex(list(1:10, NULL), index = c(1L, 1L, 2L))
    expect_identical(exp, obj@lazyIndex)

    exp <- .LazyIndex(list(1:10), index = c(1L, 1L))
    expect_identical(exp, obj0@lazyIndex)

    ## .remove_list_element
    obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj[["da0"]] <- NULL
    exp <- .LazyIndex(list(NULL), index = 1L)
    expect_identical(exp, obj@lazyIndex)

    obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj <- obj[1:10,]
    obj[["da0"]] <- NULL
    exp <- .LazyIndex(list(1:10), index = 1L)
    expect_identical(exp, obj@lazyIndex)

    ## .replace_list_element
    obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj[["da0"]] <- sample(obj[["da0"]])
    expect_identical(nullList, obj@lazyIndex)

    obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj[["da0"]] <- da0[sample(nrow(obj))]
    expect_identical(nullList, obj@lazyIndex)

    obj <- DelayedDataFrame(letters, da0 = I(da0))
    obj <- obj[10:1,]
    obj[["da0"]] <- da0[sample(10)]
    exp <- .LazyIndex(list(10:1, NULL), index=1:2)
    expect_identical(exp, obj@lazyIndex)
})

## context(".update_index")

## test_that(".update_index works", {
##     expect_identical(.LazyIndex(), .update_index(.LazyIndex(), 0, NULL))

##     lst <- .LazyIndex(list(NULL), index=1L)
##     expect_identical(lst, .update_index(lst, 1, NULL))

##     exp <- .LazyIndex(list(NULL)
## })
