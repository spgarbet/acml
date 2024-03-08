# aclme methods and utilities

# Most of this code was lifted from the lme4 package by Douglas Bates and adapted as needed 
# for working with ascertainment corrected models.


######################################################################
### Utilities for parsing the mixed model formula

### Return the pairs of expressions that separated by vertical bars
findbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("("))           return(findbars(term[[2]]))

    if (!is.call(term)) stop("term must be of class call")

    if (term[[1]] == as.name('|'))           return(term)
    if (length(term) == 2)                   return(findbars(term[[2]]))

    c(findbars(term[[2]]), findbars(term[[3]]))
}

### Return the formula omitting the pairs of expressions that are
### separated by vertical bars
nobars <- function(term)
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

### Substitute the '+' function for the '|' function
subbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
        term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

### Substitute any names from nlist in term with 1
subnms <- function(term, nlist)
{
    if (!is.language(term)) return(term)
    if (is.name(term)) {
        if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
        return(term)
    }
    stopifnot(length(term) >= 2)
    for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
    term
}

### Return the list of '/'-separated terms in an expression that
### contains slashes
slashTerms <- function(x)
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

### from a list of length 2 return recursive interaction terms
makeInteraction <- function(x)
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

### expand any slashes in the grouping factors returned by findbars
expandSlash <- function(bb)
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}