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
