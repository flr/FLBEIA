## helper function to collapse fleets.ctrl object to a data.frame
## easier to access elements but may be redundant

parse_fleets_ctrl <- function(fleets.ctrl, fleets, biols){
  flnms <- names(fleets)
  df <- do.call("rbind",lapply(flnms, function(x){
    df <- data.frame("fleet" = x,
                     "effort.model" = fleets.ctrl[[x]]$effort.model,
                     "restriction" = fleets.ctrl[[x]]$restriction,
                     "effort.restr" = fleets.ctrl[[x]]$effort.restr,
                     "stock" = unlist(lapply(names(fleets.ctrl[[x]]),function(y){
                       if("catch.model" %in% names(fleets.ctrl[[x]][[y]])){
                         y
                       } else NULL
                     })),
                     "catch.model" = unlist(lapply(names(fleets.ctrl[[x]]),function(y){
                       if("catch.model" %in% names(fleets.ctrl[[x]][[y]])){
                         fleets.ctrl[[x]][[y]]$catch.model
                       } else NULL
                     })),
                     "stock.restr" = F)
    if("stocks.restr" %in% names(fleets.ctrl[[x]])){
      df$stock.restr[df$stock %in% unlist(fleets.ctrl[[x]]$stocks.restr)] <- T
    } else {
      df$stock.restr <- T
    }
    df
  }))
  df$idx   <-  match(df$stock, names(biols))-1
  df$idx[!df$stock.restr] <- NA 
  df$fleet <- factor(df$fleet, levels = unique(df$fleet))
  df
}

## Auxiliary function to extract inputs from biols and fleets objects (stock numbers, natural mortality, lan & dis weights, selection, catchability, quota, effort shares) 
## data is organized by stock which eases the extraction of the input data when identifying a choke stock for each fleet 
## return list with two items 
##     - stock_inputs [list] dim [stocks]
##              * stock 1 [list] dim [8]
##                    - N-at-age [vector]  dim [ages]
##                    - M-at-age [vector] dim [ages]
##                    - Quota-by-fleet [vector] dim [n_fleets]
##                    - q [list] dim [fleets]
##                          ~ catchabilities [matrix] dim [metiers x ages]
##                    - cw [list] dim [fleets]
##                          ~ catch weights [matrix] dim [metiers x ages]
##                    - lw [list] dim [fleets]
##                          ~ landing weights [matrix] dim [metiers x ages]
##                    - dw [list] dim [fleets]
##                          ~ discard weights [matrix] dim [metiers x ages]
##                    - ret [list] dim [fleets]
##                          ~ retention [matrix] dim [metiers x ages]
##     - fleet_inputs [list] dim [4]
##              * eshare [vector] with effort shares [Eshare_f1_m1, Eshare_f1_m2, ..., Eshare_fn_mn]
##              * effort [vector] with effort dim [n_fleets]
##              * fl_idx [vector] fleet index 
##              * n_fl_mt [integer] number of fleet_métier combinations


get_Baranov_inputs <- function(biols,                            ## biols input to FLBEIA function
                               fleets,                           ## fleets input to FLBEIA function
                               advice,                           ## advice input to FLBEIA function
                               year = year,                      ## year input in FLBEIA simulation loop (input argument to fleets.om) 
                               unit = NULL,                      ## currently not implemented, leave NULL
                               season = NULL,                    ## currently not implemented, leave NULL
                               area = NULL,                      ## currently not implemented, leave NULL
                               iter = iter){                     ## iteration, one "input object" will be generated for each iteration
  
  if(is.null(unit)) unit_idx <- 1
  if(is.null(season)) season_idx <- 1
  if(is.null(area)) area_idx <- 1
  
  # Quota by fleet
  Qs <- advice$TAC@.Data[,year,unit_idx,season_idx,area_idx,iter]
  
  qs <- lapply(advice$quota.share, function(x){
    qs <- x@.Data[,year,unit_idx,season_idx,area_idx,iter]
    qs[qs == 0] <- NA
    names(qs) <- dimnames(x)[[1]]
    qs
  })
  
  quota_by_fleet <- lapply(names(qs), function(x){
    qs <- qs[[x]]*Qs[x]
    qs
  })
  names(quota_by_fleet) <- names(qs)
  
  # Stock numbers at the start of the year
  Ns <- lapply(biols, function(x){
    res <- x@n@.Data[,year,unit_idx,season_idx,area_idx,iter] 
    names(res) <- dimnames(x)[[1]]
    res
  })
  
  # Natural mortality
  Ms <- lapply(biols, function(x){
    res <- x@m@.Data[,year,unit_idx,season_idx,area_idx,iter] 
    names(res) <- dimnames(x)[[1]]
    res
  })
  
  # loop over stocks and add catchabilities, weights and retention from fleet object
  inputs <- vector(mode = "list", length = length(Ns))
  names(inputs) <- names(Ns)
  
  for(i in names(inputs)){
    input_names <- c("N","M","Q","q","cw","lw","dw","l.sel")
    input <- vector(mode = "list", length = length(input_names))
    names(input) <- input_names
    
    # add stock numbers, natural mortality and fleet quotas for stock i
    input$N  <- Ns[[i]]                      
    input$M  <- Ms[[i]]                      
    input$Q  <- quota_by_fleet[[i]]          
    
    # extract catchability, weights and discard selectivity from fleet objects
    flnms <- names(quota_by_fleet[[i]])
    cw.ls <- dw.ls <- lw.ls <- q.ls <- l.sel.ls <- vector(mode = "list", length = length(flnms)) 
    names(cw.ls) <- names(dw.ls) <- names(lw.ls) <- names(q.ls) <- names(l.sel.ls) <- flnms
    
    for(j in flnms){
      # métier names
      mtnms <- names(fleets[[j]]@metiers)
      
      # create matrices/vectors to store catchabilities, selectivities and weights
      cw <- lw <- dw <- q <- l.sel <- matrix(0, length(mtnms), length(Ns[[i]]), dimnames = list(mtnms, names(Ns[[i]])))
      
      
      ## catchabilities
      qsm   <- do.call("rbind",lapply(fleets[[j]]@metiers, function(x){
        if(i %in% names(x@catches)){
          x@catches[[i]]@catch.q@.Data[,year,unit_idx,season_idx,area_idx,iter]
          res <- x@catches[[i]]@catch.q@.Data[,"2023",1,1,1,1]
          names(res) <- colnames(q)
          res
        } else {
          array(0, dim = c(ncol(q)), dimnames = list("age" = colnames(q)))
        }
      }))
      
      lwm   <- do.call("rbind",lapply(fleets[[j]]@metiers, function(x){
        if(i %in% names(x@catches)){
          res <- x@catches[[i]]@landings.wt@.Data[,year,unit_idx,season_idx,area_idx,iter]
          names(res) <- colnames(q)
          res
        } else {
          array(0, dim = c(ncol(q)), dimnames = list("age" = colnames(q)))
        }
      }))
      
      dwm <- do.call("rbind",lapply(fleets[[j]]@metiers, function(x){
        if(i %in% names(x@catches)){
          res <- x@catches[[i]]@discards.wt@.Data[,year,unit_idx,season_idx,area_idx,iter]
          names(res) <- colnames(q)
          res
        } else {
          array(0, dim = c(ncol(q)), dimnames = list("age" = colnames(q)))
        }
      }))
      
      selm <- do.call("rbind",lapply(fleets[[j]]@metiers, function(x){
        if(i %in% names(x@catches)){
          res <- x@catches[[i]]@landings.sel@.Data[,year,unit_idx,season_idx,area_idx,iter]
          names(res) <- colnames(q)
          res
        } else {
          array(0, dim = c(ncol(q)), dimnames = list("age" = colnames(q)))
        }
      }))
      
      cwm <- lwm * selm + dwm * (1 - selm)      
      
      # map to correct ages
      q[,colnames(qsm)]      <- qsm
      cw[, colnames(cwm)]    <- cwm
      lw[, colnames(lwm)]    <- cwm
      dw[, colnames(dwm)]    <- cwm
      l.sel[, colnames(cwm)] <- selm
      
      cw.ls[[j]]    <- cw  
      dw.ls[[j]]    <- dw
      lw.ls[[j]]    <- lw
      q.ls[[j]]     <- q
      l.sel.ls[[j]] <- l.sel
      
    }
    input$q     <- do.call("rbind",q.ls)
    input$cw    <- do.call("rbind",cw.ls)
    input$lw    <- do.call("rbind",lw.ls)
    input$dw    <- do.call("rbind",dw.ls)
    input$l.sel <- do.call("rbind",l.sel.ls)
    
    inputs[[i]] <- input

  }
  
  ## effort
  
  effort  <- unlist(lapply(fleets, function(x) x@effort@.Data[1, year,unit_idx,season_idx,area_idx,iter]))
  n_fl_mt <- nrow(input$q)
  eshare  <- unlist(lapply(fleets, function(x){lapply(x@metiers, function(y){y@effshare@.Data[1,year,unit_idx,season_idx,area_idx,iter]})}))
  fl_idx  <- unlist(lapply(1:length(fleets), function(x) rep((x-1), length(fleets[[x]]@metiers))))
  
  return(list("stock_inputs" = inputs, "fleet_inputs" = list("eshare" = eshare, "effort" = effort, "fl_idx" = fl_idx, "n_fl_mt" = n_fl_mt)))
}


find_effort_choke_Baranov <- function(idx, baranov_inputs, fleets.ctrl.ls, estart = NULL){
  
  #### STEP 1: Identify potential choking stocks for each fleet
  
  # does a fleet has a non NA quota
  qs               <- do.call('cbind',lapply(baranov_inputs$stock_inputs, function(x)x$Q))
  indices          <- apply(qs,1,function(x)(which(!is.na(x))-1))
  
  # if "stock.restr" are present, select only stocks within the "stock.restr" vector for each fleet
  indices <- lapply(names(indices), function(x){
    intersect(indices[[x]],fleets.ctrl.ls[[x]]$idx)
  })
  
  
  # modify initial idx if a stock is selected for a fleet that cannot choke (based on the indices list created above) 
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  idx <- unlist(lapply(1:length(idx), function(x){
    if(!(idx[x]%in%indices[[x]])){
      idx[x] <- resample(indices[[x]],1)
    }
    idx[x]
  }))
  
  ## modify the index of constant effort fleets to -1
  cstE_idx         <- unlist(lapply(fleets.ctrl.ls,function(x) unique(x$effort.model))) == "fixedEffort"
  idx[cstE_idx]    <- -1
  
  #### STEP 2: PREPARE INPUTS FOR SOLVER
  
  ## 2.1 identify initial vector of effort values (exclude "fixedEffort" fleets) 
  if(is.null(estart)) estart <- baranov_inputs$fleet_inputs$effort[idx >= 0] 
  if(length(estart) != sum(idx>=0)) estart <- baranov_inputs$fleet_inputs$effort[idx >= 0] 
  
  ## 2.2 select catch of landing weights for the inputs based on the restriction, modify catchabilities in case of landings
  restricition   <- unlist(lapply(fleets.ctrl.ls,function(x)unique(x$restriction)))
  nleqslv_inputs <- baranov_inputs
  if(all(restricition == "catch")){
    nleqslv_inputs$stock_inputs <- lapply(nleqslv_inputs$stock_inputs, function(x){
      x$wt    <- x$cw
      x$q.sel <- x$q
      x$cw    <- x$lw <- x$dw <- x$l.sel <- NULL
      x
    }) 
  } else {
    nleqslv_inputs$stock_inputs <- lapply(nleqslv_inputs$stock_inputs, function(x){
      x$wt    <- x$cw
      x$q.sel <- x$q
      for(y in 1:length(restricition)){
        if(restriction[y] == "landings"){
          f_idx <- length(unlist(nleqslv_inputs$fleet_inputs$effortshare[0:(y-1)])) + (1:length(nleqslv_inputs$fleet_inputs$effortshare[[y]]))
          x$wt[f_idx]    <- x$lw[f_idx]
          x$q.sel[f_idx] <- x$q[f_idx] * x$l.sel[f_idx]
        } 
      }
      x$cw <- x$lw <- x$dw <- x$l.sel <- NULL
      x
    })
  }
  
  
  
  # solve system of nonlinear eqauations
  res <- nleqslv(estart, Q2baranov, inputs = nleqslv_inputs, idx = idx, method = "Newton", global = "qline")
  
  
  #### STEP 3: IDENTIFY SMFB CONDITIONS
  SMFB_fleets <- unlist(lapply(fleets.ctrl.ls,function(x)unique(x$effort.model))) == "SMFB"
  min_fleets  <- unlist(lapply(fleets.ctrl.ls[SMFB_fleets],function(x)(unique(x$effort.restr) == "min")))
  max_fleets  <- unlist(lapply(fleets.ctrl.ls[SMFB_fleets],function(x)(unique(x$effort.restr) == "max")))
  
  
  ## compute catch for all fleets, except fixedEffort
  
  catch <- do.call('rbind',lapply(0:(length(baranov_inputs$stock_inputs)-1), function(x) {
    idx_stock_x <- rep(x,length(idx))
    idx_stock_x[idx < 0] <- -1
    Q2baranov(res$x, nleqslv_inputs, idx_stock_x)
  }))
  quota <- do.call("rbind",lapply(nleqslv_inputs$stock_inputs,"[[","Q"))[,SMFB_fleets]
  dimnames(catch) <- dimnames(quota)
  
  ## identify stop conditions
  ## compare catch to quota of SMFB fleets
  ## only consider restricting stocks
  any_quota_overshot <- any_quota_undershot <- F
  
  catch_check   <- catch
  catch_check[] <- NA
  flnms         <- names(fleets.ctrl.ls)
  for(i in 1:length(fleets.ctrl.ls)){  
    if(SMFB_fleets[i]){
      catch_check[fleets.ctrl.ls[[i]]$idx[!is.na(fleets.ctrl.ls[[i]]$idx)] + 1,flnms[i]] <- catch[fleets.ctrl.ls[[i]]$idx[!is.na(fleets.ctrl.ls[[i]]$idx)] + 1,flnms[i]]
    }
  }
  
  if(any(min_fleets)) any_quota_overshot      <- any(apply(catch_check > 1e-6, 2, any, na.rm = T)[min_fleets])
  if(any(max_fleets)) any_quota_undershot     <- any(apply(catch_check < -1e-6, 2, any, na.rm = T)[max_fleets])
  
  
  # update solution if necessary
  ctr <- 0
  while(any_quota_overshot | any_quota_undershot){
    ctr <- ctr + 1
    print(ctr)
    
    ## min fleets
    ## replace current stock tested for choking (if any catch > quota) by stock with highest overshoot
    modif <- colSums(catch_check > 1e-6, na.rm = T)>0 & min_fleets
    idx[idx>=0][modif] <- apply((catch_check + quota)/quota, 2, function(x)(which.max(x)-1))[modif]
    
    ## max fleets
    ## replace current stock tested for choking (if any catch < quota) by stock with highest undershoot
    modif <- colSums(catch_check < -1e-6, na.rm = T)>0 & max_fleets
    idx[idx>=0][modif] <- apply((catch_check + quota)/quota, 2, function(x)(which.min(x)-1))[modif]
    
    
    # use previous effort as start value
    estart <- res$x
    res    <- nleqslv(estart, Q2baranov, control=list(btol=.01), inputs = nleqslv_inputs, idx = idx, method = "Newton", global = "qline")
    
    # compute quota minus catch to compute stop conditions
    catch <- do.call('rbind',lapply(0:(length(stocks)-1), function(x) {
      idx_stock_x <- rep(x,length(idx))
      idx_stock_x[idx < 0] <- -1
      Q2baranov(res$x, nleqslv_inputs, idx_stock_x)
    }))
    dimnames(catch) <- dimnames(quota)
    
    ## only consider restricting stocks
    catch_check   <- catch
    catch_check[] <- NA
    flnms         <- names(fleets.ctrl.ls)
    for(i in 1:length(fleets.ctrl.ls)){  
      if(SMFB_fleets[i]){
        catch_check[fleets.ctrl.ls[[i]]$idx[!is.na(fleets.ctrl.ls[[i]]$idx)] + 1,flnms[i]] <- catch[fleets.ctrl.ls[[i]]$idx[!is.na(fleets.ctrl.ls[[i]]$idx)] + 1,flnms[i]]
      }
    }
    if(any(min_fleets)) any_quota_overshot      <- any(apply(catch_check > 1e-6, 2, any, na.rm = T)[min_fleets])
    if(any(max_fleets)) any_quota_undershot     <- any(apply(catch_check < -1e-6, 2, any, na.rm = T)[max_fleets])
  }
  cat(paste0("nleqslv termination code: ", res$termcd))
  cat("\n")
  
  # return effort by fleet
  effort           <- nleqslv_inputs$fleet_inputs$effort  # copy fixedEffort for all fleets
  effort[idx >= 0] <- res$x                               # replace fixedEffort with new effort for SMFB fleets
  return(list("effort" = effort, "choke" = idx, "catch" = catch + quota, "quota_upt" = round(100*(catch + quota)/quota, 2)))
}



## call nleqslv solver
solve_Baranov <- function(biols, fleets, fleets.ctrl.df, advice, year, baranov_inputs_ls){
  
  fleets.ctrl.ls   <- split(fleets.ctrl.df, fleets.ctrl.df$fleet)
  
  lapply(baranov_inputs_ls, function(x){
    if(any(fleets.ctrl.df$effort.model == "SMFB")){
      idx              <- unlist(lapply(fleets.ctrl.ls, function(x)x$idx[!is.na(x$idx)][1]))
      idx[unlist(lapply(fleets.ctrl.ls, function(x)x$effort.model[1])) == "fixedEffort"] <- -1
      ## initial vector of fishing effort by fleet
      estart           <- rep(0.01,sum(idx>=0))
      out <- find_effort_choke_Baranov(idx = idx, baranov_inputs = x, fleets.ctrl.ls = fleets.ctrl.ls, estart = estart)
      
    }
    if(all(fleets.ctrl.df$effort.model == "fixedEffort")){
      out <- list()
      out$effort <- x$fleet_inputs$effort
    }
    out
  })
}

Baranov_simulatneous <- function(biols = biols, fleets = fleets, BDs = BDs, advice = advice,
                                 year = year, season = season, biols.ctrl=biols.ctrl, fleets.ctrl = fleets.ctrl, covars = covars, 
                                 assess.ctrl=assess.ctrl, advice.ctrl = advice.ctrl){
  nit <- dim(biols[[1]]@n)[6]
  
  fleets.ctrl.df    <-  parse_fleets_ctrl(fleets.ctrl, fleets, biols)
  baranov_inputs_ls <- lapply(1:nit, function(i) get_Baranov_inputs(biols, fleets, advice, year = year, iter = i))
  res               <- solve_Baranov(biols, fleets, fleets.ctrl.df, advice, year, baranov_inputs_ls)
  tmp_effort        <- do.call("rbind", lapply(res,"[[","effort"))
  fleets            <- update_catch_effort_FLFleets(tmp_effort, fleets, baranov_inputs_ls, (year-1))
  list(fleets = fleets)
}