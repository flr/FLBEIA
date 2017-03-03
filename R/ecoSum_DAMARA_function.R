#------------------------------------------------------------------------------#
#    Additional summary function for DAMARA taking account of
#    new covars for fleet economics
#   
# Paul Dolder / Simon Mardle
# Created: 14/10/2015
# Changed: 
#------------------------------------------------------------------------------#

# Requires Covar inputs
# FuelCost - per unit of effort
# CapitalCost - per vessel
# Salaries - blank, is calculated
# NumbVessels -
# MaxDays    - capacity (in unit of effort)
# MaxDaysPerVessel - capacity per vessel (in unit of effort)
# EmploymentPerVessel - number crew per vessel
# VariableCost  - per unit of effort
# FixedCost - per vessel
# InvestShare - proportion
# OtherRevenueFocusArea - revenue from other species, calculated in projections as per unit of effort multiplier
# OtherRevenueElsewhere - revenue for fleets outside of model bounds - calculated in projections as fixed per vessel
# DepreciationCost - per vessel
# CapitalValue - per vessel
# KWdaysPerVessel - effort per vessel
# AvgKwPerVessel - average kw power of a vessel in fleet
# KWDaysFocusArea - effort per vessel spent on focus area

#------------------------------------------------------------------------------#
# ecoSum data.frame[year, quarter, stock, fleet, iter, ||,|| 
#        profits, capacity, costs, discards, effort, landings,
#        DAS_FocusArea, DAS_Elsewhere, revenueFocusArea, revenueElsewhere,
#       totalRevenue, crewCosts, fuelCosts, variableCosts,fixedCosts,
#       depreciationCosts, investmentCosts,GCF,GVA,netProfit,BER,employment,
#       numberVessels] 
#------------------------------------------------------------------------------#
#' @rdname summary_flbeia
ecoSum_damara <- function (fleets, flnms = "all", years, covars = NULL)
{
    if (flnms[1] == "all")
        flnms <- names(fleets)
    Dim <- dim(fleets[[1]]@effort[, years, ])[c(2, 4, 6)]
    Dimnm <- dimnames(fleets[[1]]@effort[, years, ])
    n <- prod(Dim) * length(flnms)
    res <- data.frame(year = rep(years, prod(Dim[2:3]) * length(flnms)),
        quarter = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3] *
            length(flnms)), fleet = rep(flnms, each = prod(Dim)),
        iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(flnms)),
        capacity = numeric(n), costs = numeric(n), effort = numeric(n),
        profits = numeric(n),               ## Adaptation from here to bring in additional covars
        DAS_FocusArea=numeric(n),DAS_Elsewhere=numeric(n),
        revenueFocusArea=numeric(n),revenueElsewhere=numeric(n),
        totalRevenue=numeric(n),crewCosts=numeric(n),fuelCosts=numeric(n),
        variableCosts=numeric(n),fixedCosts=numeric(n),
        depreciationCosts=numeric(n),investmentCosts=numeric(n),
        GCF=numeric(n),GVA=numeric(n),netProfit=numeric(n),BER=numeric(n),
        employment=numeric(n),numberVessels=numeric(n))
    k <- 1
    for (f in flnms) {
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        res[k:(k + prod(Dim) - 1), "capacity"] <- c(fl@capacity[,
            years, ])
        res[k:(k + prod(Dim) - 1), "effort"] <- c(fl@effort[,
            years, ])
        if(!is.null(covars)) res[k:(k + prod(Dim) - 1), "costs"] <- c(costs_flbeia(fl, covars, f)[,years, ])
        res[k:(k + prod(Dim) - 1), "profits"] <- c(revenue_flbeia(fl)[,
            years, ]) - res[k:(k + prod(Dim) - 1), "costs"]
        ## Additional DAMARA outputs    
        if(!is.null(covars)) {
        ## fleet based costs (i.e. fixed and capital)
        res[k:(k + prod(Dim) - 1), "fixedCosts"] <- res[k:(k + prod(Dim) - 1), "fixedCosts"] + c(covars[["FixedCost"]][f,years,] * covars[["NumbVessels"]][f, years,])
        res[k:(k + prod(Dim) - 1), "depreciationCosts"] <- res[k:(k + prod(Dim) - 1), "depreciationCosts"] + c(covars[["DepreciationCost"]][f, years,] * covars[["NumbVessels"]][f, years,])
        res[k:(k + prod(Dim) - 1), "investmentCosts"] <- res[k:(k + prod(Dim) - 1), "investmentCosts"] + c(covars[["InvestShare"]][f, years,] * covars[["CapitalValue"]][f, years,] * covars[["NumbVessels"]][f, years,])
        ##metier based costs (i.e. fuel and other variable costs)
        for (mt in mts) {   ##CHECK ALL EFFORT IS ACCOUNTED FOR!
          res[k:(k + prod(Dim) - 1), "fuelCosts"] <- res[k:(k + prod(Dim) - 1), "fuelCosts"] + c(covars[["FuelCost"]][f, years,] * fl@effort[,years,] * fl@metiers[[mt]]@effshare[,years,])
          res[k:(k + prod(Dim) - 1), "variableCosts"] <- res[k:(k + prod(Dim) - 1), "variableCosts"] + c(covars[["VariableCost"]][f, years,] * fl@effort[,years,] * fl@metiers[[mt]]@effshare[,years,])
          #calc DAS elsewhere???
          res[k:(k + prod(Dim) - 1), "DAS_FocusArea"] <- res[k:(k + prod(Dim) - 1), "DAS_FocusArea"] + c(fl@effort[,years,] * fl@metiers[[mt]]@effshare[,years,]) / c(covars[["NumbVessels"]][f,years,]*covars[["AvgKwPerVessel"]][f,years,])
          ##revenues
          m <- fl@metiers[[mt]]
          sts <- catchNames(fl)
          for (st in sts) {
             if (!(st %in% catchNames(m))) 
                next
             dat <- m@catches[[st]]
             res[k:(k + prod(Dim) - 1), "revenueFocusArea"] <- res[k:(k + prod(Dim) - 1), "revenueFocusArea"] + 
                                c(apply(dat@landings.n[,years,] * dat@landings.wt[,years,] * dat@price[,years,], c(2, 4, 6), sum, na.rm = T))  
          }
        }
        ## Need to add in additional revenue from outside 7B-K and other species - check FG and BK
        res[k:(k + prod(Dim) - 1), "revenueFocusArea"] <- res[k:(k + prod(Dim) - 1), "revenueFocusArea"] + c(covars[["OtherRevenueFocusArea"]][f, years,] * fl@effort[,years,] * fl@metiers[[mt]]@effshare[,years,])
        res[k:(k + prod(Dim) - 1), "revenueElsewhere"] <- res[k:(k + prod(Dim) - 1), "revenueElsewhere"] + c(covars[["NumbVessels"]][f, years,] * covars[["OtherRevenueElsewhere"]][f, years,])
        res[k:(k + prod(Dim) - 1), "totalRevenue"] <- (res[k:(k + prod(Dim) - 1), "revenueFocusArea"] + res[k:(k + prod(Dim) - 1), "revenueElsewhere"])
        ##crewCosts
        res[k:(k + prod(Dim) - 1), "crewCosts"] <- res[k:(k + prod(Dim) - 1), "totalRevenue"] * c(fl@crewshare[,years,])
        ##profits and BER
        res[k:(k + prod(Dim) - 1), "GCF"] <- res[k:(k + prod(Dim) - 1), "totalRevenue"] - res[k:(k + prod(Dim) - 1), "fixedCosts"] - res[k:(k + prod(Dim) - 1), "fuelCosts"] - res[k:(k + prod(Dim) - 1), "variableCosts"] - res[k:(k + prod(Dim) - 1), "crewCosts"]
        res[k:(k + prod(Dim) - 1), "GVA"] <- res[k:(k + prod(Dim) - 1), "GCF"] + res[k:(k + prod(Dim) - 1), "crewCosts"]
        res[k:(k + prod(Dim) - 1), "netProfit"] <- res[k:(k + prod(Dim) - 1), "GCF"] - res[k:(k + prod(Dim) - 1), "depreciationCosts"] 
        res[k:(k + prod(Dim) - 1), "BER"] <- (res[k:(k + prod(Dim) - 1), "crewCosts"] + res[k:(k + prod(Dim) - 1), "fixedCosts"] + res[k:(k + prod(Dim) - 1), "depreciationCosts"])/
                                             (1-(res[k:(k + prod(Dim) - 1), "fuelCosts"]/ res[k:(k + prod(Dim) - 1), "totalRevenue"]) - (res[k:(k + prod(Dim) - 1), "variableCosts"])/res[k:(k + prod(Dim) - 1), "totalRevenue"])
        res[k:(k + prod(Dim) - 1), "employment"] <- res[k:(k + prod(Dim) - 1), "employment"] + c(covars[["EmploymentPerVessel"]][f,years,] * covars[["NumbVessels"]][f,years,])
        res[k:(k + prod(Dim) - 1),"numberVessels"] <-res[k:(k + prod(Dim) - 1), "numberVessels"] + c(covars[["NumbVessels"]][f,years,])
        }
        k <- k + prod(Dim)
    }
    return(res)
}
