###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        xtabs2
# NOTE #1:      Return a contingency table as an array, similar to xtabs.
#           xtabs converts NA-s in 0, but xtabs2 does not touch NA values.
#           The class of the output in xtabs is "xtabs" "table" 
#           and in the function xtabs2 is an array.
###############################################################################

  
  xtabs2 <- function (formula = ~., data = parent.frame(), subset, sparse = FALSE, 
                      na.action, exclude = c(NA, NaN), drop.unused.levels = FALSE) 
  {
    if (missing(formula) && missing(data)) 
      stop("must supply either 'formula' or 'data'")
    if (!missing(formula)) {
      formula <- as.formula(formula)
      if (!inherits(formula, "formula")) 
        stop("'formula' missing or incorrect")
    }
    if (any(attr(terms(formula, data = data), "order") > 1)) 
      stop("interactions are not allowed")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
      m$data <- as.data.frame(data)
    m$... <- m$exclude <- m$drop.unused.levels <- m$sparse <- NULL
    m[[1L]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    if (length(formula) == 2L) {
      by <- mf
      y <- NULL
    }
    else {
      i <- attr(attr(mf, "terms"), "response")
      by <- mf[-i]
      y <- mf[[i]]
    }
    by <- lapply(by, function(u) {
      if (!is.factor(u)) 
        u <- factor(u, exclude = exclude)
      u[, drop = drop.unused.levels]
    })
    if (!sparse) {
      x <- if (is.null(y)) 
        do.call("table", by)
      else if (NCOL(y) == 1L) 
        tapply(y, by, sum)
      else {
        z <- lapply(as.data.frame(y), tapply, by, sum)
        array(unlist(z), dim = c(dim(z[[1L]]), length(z)), 
              dimnames = c(dimnames(z[[1L]]), list(names(z))))
      }
      x
    }
    else {
      if (length(by) != 2L) 
        stop("xtabs(*, sparse=TRUE) applies only to two-way tables")
      if (is.null(tryCatch(loadNamespace("Matrix"), error = function(e) NULL))) 
        stop("xtabs(*, sparse=TRUE) needs package \"Matrix\" correctly installed")
      if (length(i.ex <- unique(unlist(lapply(by, function(f) which(is.na(f))))))) 
        by <- lapply(by, `[`, -i.ex)
      rows <- by[[1L]]
      cols <- by[[2L]]
      rl <- levels(rows)
      cl <- levels(cols)
      if (is.null(y)) 
        y <- rep.int(1, length(rows))
      as(new("dgTMatrix", i = as.integer(rows) - 1L, j = as.integer(cols) - 
               1L, x = as.double(y), Dim = c(length(rl), length(cl)), 
             Dimnames = list(rl, cl)), "CsparseMatrix")
    }
  }