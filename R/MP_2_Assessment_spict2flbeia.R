


spict2flbeia <- function(stock, indices, control=NULL,covars=covars){
  # Hack input
  
  # check if spict results slot already exist
  # if they are missing, add them
  slot_names <- c(paste0("spict_", c("F", "B", "Fmsy", "Bmsy",
                                     "FFmsy", "BBmsy")))
  st <- name(stock)
  results <- list()
  years <- dimnames(stock@stock)$year
  res_template <- stock@stock[,years] 
  res_template[] <- NA
  for(j in 1:length(slot_names)){
    results[[st]][[j]]<- res_template
    print(j)
  }
  names(results[[st]]) <- slot_names
  
  if(any(slot_names %in% names(covars[[st]]))){
    for(j in 1:2){
      first.yr <- as.numeric(years[1])
      last.yr <- as.numeric(tail(years,n=1))
      covars[[st]][[j]] <- window(covars[[st]][[j]],first.yr,last.yr)}
    }else{
       # create template FLQuant for storing assessment results
       # add slots
       covars <- list()
       covars[[st]] <- list()
       for(j in 1:2){
         covars[[st]][[j]]<- res_template
       }
       names(covars[[st]]) <- slot_names[c(3,4)]
 }
  
  # Iterations
  
  for (i in 1:dim(catch(stock))[6]){
    print(i)
    ip_list <- vector("list")
    ip_list$obsC <- c(iter(catch(stock),i))
    ip_list$timeC <- as.numeric(dimnames(catch(stock))$year)
    ip_list$timeI <- lapply(indices, function(x) as.numeric(dimnames(index(x))$year))
    ip_list$obsI <- lapply(indices, function(x) c(iter(index(x),i)))
    ip_list$dteuler <- 1/8
    ip_list$getReportCovariance <- FALSE
    # Fit
    out <- fit.spict(ip_list)
    #sumspict.states(out)
    
    ### extract output
    # get desired years
    years_wanted <- dimnames(catch(stock))$year
    years_available <- row.names(get.par("logB", out, exp = TRUE))
    years <- intersect(years_wanted, years_available)
    
    # extract values
    pars_des <- c(paste0("log", c("F", "B", "Fmsy", "Bmsy", "FFmsy", "BBmsy")))
    # extract times series
    for(j in seq_along(pars_des[c(1:2, 5:6)])){
      res <- get.par(pars_des[c(1:2, 5:6)][j], out, exp = TRUE)[, "est"]
      # subset to desired years
      res <- res[names(res) %in% years]
      # insert into stock
      iter(results[[st]][[slot_names[c(1:2, 5:6)[j]]]][, years], i)[] <- res
    }
    # extract reference points
    for(j in seq_along(pars_des[c(3:4)])){
      res <- get.par(pars_des[c(3:4)][j], out, exp = TRUE)[, "est"]
      # insert into stock
      iter(results[[st]][[slot_names[c(3:4)[j]]]][, tail(dimnames(stock@stock)$year,1)], i)[] <- res
    }
    
  }
  
  stock@stock[,years] <- results[[st]][["spict_B"]] #explotable biomass
  stock@stock.n[,years] <- stock@stock/stock@stock.wt
  stock@harvest[,years] <- results[[st]][["spict_F"]] #explotable biomass
  covars[[st]]$spict_Bmsy[,tail(years,n=1)] <- results[[st]][["spict_Bmsy"]][,tail(years,n=1)]
  covars[[st]]$spict_Fmsy[,tail(years,n=1)] <- results[[st]][["spict_Fmsy"]][,tail(years,n=1)]
  # Returns what?
  # Depends on HCR
  # Update stock with biomass estimates
  rm(list=setdiff(ls(), c("stock","covars")))
  return(list(stock = stock,covars=covars))
}

