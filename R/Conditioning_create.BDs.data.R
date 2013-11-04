
 create.BDs.data <- function (path)
{
    nmstks <- unique(sub(".*?^(.*?)_bd.model*", "\\1", ls(pattern = "_bd.model",
        envir = sys.frame(which = 0))))
    n.stk.BD <- length(nmstks)
    if (n.stk.BD != 0) {
        list.stks.flq <- create.list.stks.flq()
        ny <- length(first.yr:last.yr)
        hist.yrs <- as.character(first.yr:(proj.yr - 1))
        proj.yrs <- as.character(proj.yr:last.yr)
        for (i in 1:n.stk.BD) {
            nmstk <- nmstks[i]
            cat("=============", nmstk, "BD", "=============\n")
            stk.index <- which(names(list.stks.flq) == nmstks)
            flq.stk <- list.stks.flq[[stk.index]][, , 1]
            stk.model <- get(paste(nmstk, "_bd.model", sep = ""))
            stk.unit <- get(paste(nmstk, ".unit", sep = ""))
            stk.biomass <- get(paste(nmstk, "_biomass.flq", sep = ""))
            stk.catch <- get(paste(nmstk, "_catch.flq", sep = ""))
            stk.param.n <- get(paste(nmstk, "_param.n", sep = ""))
            stk.params.name <- get(paste(nmstk, "_params.name",
                sep = ""))
            stk.params <- get(paste(nmstk, "_params.array", sep = ""))
            stk.range.min <- get(paste(nmstk, "_range.min", sep = ""))
            stk.range.max <- get(paste(nmstk, "_range.max", sep = ""))
            stk.range.plusgroup <- get(paste(nmstk, "_range.plusgroup",
                sep = ""))
            stk.range.minyear <- get(paste(nmstk, "_range.minyear",
                sep = ""))
            stk.uncertainty <- mget(paste(nmstk, "_uncertainty.flq",
                sep = ""), envir = as.environment(-1), ifnotfound = NA,
                inherits = TRUE)[[1]]
            params <- array(dim = c(stk.param.n, ny, ns, ni),
                dimnames = list(param = ac(1:stk.param.n), year = ac(first.yr:last.yr),
                  season = ac(1:ns), iter = 1:ni))
            stk.bd <- FLBDsim(name = nmstk, model = stk.model,
                biomass = flq.stk, catch = flq.stk, uncertainty = flq.stk,
                params = params)
            dimnames(stk.bd@params)$param <- stk.params.name
            stk.bd@range[["min"]] <- stk.range.min
            stk.bd@range[["max"]] <- stk.range.max
            stk.bd@range[["plusgroup"]] <- stk.range.plusgroup
            stk.bd@range[["minyear"]] <- stk.range.minyear
            stk.bd@range[["maxyear"]] <- proj.yr - 1
            if (!all(is.na(stk.uncertainty))) {
                log.dim <- equal.flq.Dimnames(lflq = list(stk.uncertainty,
                  stk.bd@uncertainty), 2)
                if (!log.dim)
                  stop("BD uncertainty dimension names \n")
                if (!(any(dim(stk.uncertainty)[3] == c(1, stk.unit))))
                  stop("in uncertainty number of stock units 1 or stk.unit")
                if (!(any(dim(stk.uncertainty)[4] == c(1, ns))))
                  stop("in uncertainty number of seasons 1 or ns")
                if (!(any(dim(stk.uncertainty)[6] == c(1, ni))))
                  stop("in uncertainty number of iterations 1 or ni")
            }
            else {
                stk.uncertainty = 1
                cat("BD uncertainty = 1 \n")
            }
            if (!all(is.na(stk.biomass))) {
                log.dim <- equal.flq.Dimnames(lflq = list(stk.biomass,
                  stk.bd@biomass[, hist.yrs]), 2)
                if (!log.dim)
                  stop("in BD biomass dimension names \n")
                if (!(any(dim(stk.biomass)[3] == c(1, stk.unit))))
                  stop("in biomass number of stock units 1 or stk.unit")
                if (!(any(dim(stk.biomass)[4] == c(1, ns))))
                  stop("in biomass number of seasons 1 or ns")
                if (!(any(dim(stk.biomass)[6] == c(1, ni))))
                  stop("in biomass number of iterations 1 or ni")
            }
            else {
                cat("BD biomass values all NA-s \n")
            }
            if (!all(is.na(stk.catch))) {
                log.dim <- equal.flq.Dimnames(lflq = list(stk.catch,
                  stk.bd@catch[, hist.yrs]), 2)
                if (!log.dim)
                  stop("in BD catch dimension names \n")
                if (!(any(dim(stk.catch)[3] == c(1, stk.unit))))
                  stop("in catch number of stock units 1 or stk.unit")
                if (!(any(dim(stk.catch)[4] == c(1, ns))))
                  stop("in catch number of seasons 1 or ns")
                if (!(any(dim(stk.catch)[6] == c(1, ni))))
                  stop("in catch number of iterations 1 or ni")
            }
            else {
                cat("BD catch values all NA-s \n")
            }
            log.dim <- equal.flq.Dimnames(lflq = list(stk.params,
                stk.bd@params), 1:4)
            if (!log.dim)
                stop("BD parameters dimension names \n")
            stk.bd@biomass[, hist.yrs] <- stk.biomass
            stk.bd@catch[, hist.yrs] <- stk.catch
            stk.bd@params <- stk.params
            dimnames(stk.bd@params)$param <- stk.params.name
            stk.bd@uncertainty[] <- stk.uncertainty
            if (!any(is.na(stk.params[, proj.yrs, , ]))) {
                if (!all(dim(stk.params) == dim(stk.bd@params))) {
                  stop("in BD parameters dimension names \n")
                }
            }
            else {
                stop("BD parameters all NA-s \n")
            }
            if (any(is.na(stk.bd@uncertainty[, proj.yrs]))) {
                stop("Na values in uncertainty in the projection years")
            }
            assign(paste(nmstk, ".bd", sep = ""), stk.bd)
        }
        stk.bds <- paste(nmstks, ".bd", sep = "")
        BDs <- sapply(stk.bds, get, envir = sys.frame(sys.parent(-1)),
            simplify = FALSE)
        names(BDs) <- nmstks
    }
    return(BDs)
}
