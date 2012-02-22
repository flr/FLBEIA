#-------------------------------------------------------------------------------
#                           SelAtAge
#
# - _SelAtAge_ is a function to estimate the historical selectivity at age.
#   The function has been built from a template created to calculate NHKE's 
#   selectivy at age. 
#
#   To calculate de selectivity we used the formula:
#
#               Cafm = (Cfm/sum(Sifm*Bi, i=a0,...,a+))*Safm*Ba
#
#   [a = age, f = fleet, m = metier, i subscript for age]  
#   Safm = Selectivity at age 'a' for fleet 'f' and metier 'm'.
#   Cafm = catch (in weight) at age 'a' for fleet 'f' and metier 'm'.
#    Cfm = Total catch for fleet 'f' and metier 'm'.
#     Ba = Biomass (in weight) at age.
#
# See the latex doc to see the derivation of the formula.
#
# The equation above is nonlinear and  we cannot find an analytical expression
# for the selectivity.  
# Rewriting the equation above for each 'a' we have de following equation:
#   sum(Sifm*Bi, i=a0,...,a+) - (Cfm/Cafm)*Ba*Safm  
# Thus we have a system of linear equations being Sa0fm,...,Sa+fm the unknown variables.
# The problem is that the equations in the system are not independent. 
# We have to remove one equation to get a system of independent equations,
# Furthermore we have to add a constraint (an equation) to be able to solve the 
# system. 
# As we have to do this in a mechanic way, and we want to be sure that we are 
# doing it OK, in its year and season, we solve the system 3 times, removing a 
# different equation each time. And we use the same contraint (equation):
# sum(Sifm, i=a0,...a+) = 1.
# Then we compare the values of Sifm and we check if they are the same.
#
# 'restriction' argument:
#   if restriction = 1 => sum(Sa) = 1.
#   if restriction = 2 => max(Sa) = 1.
# 
# Created: 20/09/2010 10:26:30
# Dorleta Garcia.
# Last change: 2011-03-18 09:46:35 (Sonia Sánchez)
#-------------------------------------------------------------------------------


SelAtAge <- function(biols, fleets, flnm = 1, mtnm = 1, stnm = 1, years, iter = NULL, restriction = 1, ntrials = 3){

    if(class(biols) == 'FLBiols')  biol <- biols[[stnm]]
    else biol <- biols
    
    if(class(fleets) == 'FLFleetsExt') flt <- fleets[[flnm]]
    else flt <- fleets
     
    mtr <- flt@metiers[[mtnm]]
    ctc <- mtr@catches[[stnm]]
    
    ns <- dim(biol@n)[4]
    if(is.null(iter)) it <- 1:dim(biol@n)[6]  else it <- iter
    
    l <- ntrials
    
    #-------------------------------------------------------------------------------
    #  Transform the data into arrays in order to have natural ages instead of
    #  ages and units
    #------------------------------------------------------------------------------

    Bay <- unit2age(biol@wt*biol@n)
    Cay <- unit2age((catch.n(ctc)*catch.wt(ctc)))   #[na*nu,ny,ns,it]
    dimnames(Bay)[[2]] <- dimnames(biol@wt)[[2]]
    dimnames(Cay)[[2]] <- dimnames(catch.n(ctc))[[2]]
    
    # Total catch.
    Cty <- apply(Cay,2:4,sum)
    
    Say0 <- Sayk <- Bay; Say0[,,,] <- Sayk[,,,] <- 0
    
    Say <- vector('list', l)
    for(k in 1:l) Say[[k]] <- Sayk

    for(j in it){
        for(y in years){
            for(s in 1:ns){
                ok.flag <- NA
        
                pos <- unique(c(which(Bay[,y,s,j] == 0), which(Cay[,y,s,j] == 0)))
            
                # If zero catches -> Sa[,y,s,] == 0
                if(length(pos) == dim(Bay)[1]) {
                    Say0[,y,s,] <- 0
                    next
                }
            
                if(length(pos) == 0) pos <-  dim(Bay)[1] + 100
                
                Ba  <- Bay[-pos,y,s,j]
                Ca  <- Cay[-pos,y,s,j]
                Ct  <- Cty[y,s,j] 
        
                b. <-  rep(0,length(Ba)) #
                A. <- diag(0,length(Ba))
        
                for(i in 1:dim(A.)[1]){
                    A.[i,i]  <- (Ba[i] - Ct*Ba[i]/Ca[i])
                    A.[i,-i] <- Ba[-i] 
                }
        
                n <- sample(1:(dim(A.)[1]),l)
            
                for (k in 1:l) {
            
                    A <- A.; b <- b.
                    A[n[k],] <- 1
                    b[n[k]] <- 1
              
                    res <- try(solve(A,b))
              
                    if(class(res) == 'try-error'){
                        cat('-------------------------------------------------------\n')
                        cat('Error in year, ', y, ', season, ', s, ' and equation', k, ' \n')
                        cat('-------------------------------------------------------\n')
                    } else{
                        Say[[k]][-pos,y,s,j] <- res
                        ok.flag <- k
    #                    print(k)
                    }
        
                    if(pos[1] < 100)
                        for (k in 1:l) Say[[k]][pos,y,s,j] <- 0 
        

                }
                
                if(is.na(ok.flag)) stop(paste('The ',l,'trial(s) to find the selectivity have failed!!!!!')) 
                else {    
                    for (k in l:1){ 
                        if(ok.flag == k) { 
                        Say0[,y,s,j] <- Say[[k]][,y,s,j]
                            
                        }
                    }
                }
            }
        } 
    }
    
    if(l > 1){ 
        dif <- c()
        for (k in 1:(l-1)) for (r in (k+1):l) dif <- c( dif, range(Say[[k]]/Say[[l]],na.rm=T)) 
    
        dif <- range(dif)
        if(dif> 0.9 & dif < 1.1){
            cat('*******************************************************\n')
            cat(' THE PROCESS SEEMS TO HAVE BEEN SUCCESSFUL \n')
            cat('*******************************************************\n')
        }
        else{
            cat('*******************************************************\n')
            cat(' IT SEEMS THERE HAVE BEEN A PROBLEM IN THE PROCESS!!!!!\n')
            cat('*******************************************************\n')        
        }
    }
    
    if(restriction == 2) # max(Say0) = 1.
        Say0 <- sweep(Say0,2:4, apply(Say0,2:4,max), "/")
    
    # They are equal. Say = Selectividad total.
    Say <- age2unit(Say0, biol@wt)

    # Then, the selectivity should be discomposed in discards selectivity and 
    # landings selectivity.
    # dis.sel = Say*(Dis.n*Dis.wt)/(Catch.n*Catch.wt)   (dis.wt=lan.wt=catch.wt => lo quitamos)
    # lan.sel = Say*(lan.n*Dis.wt)/(Catch.n*Catch.wt)

    
    return(Say)
}
      
                                                                    