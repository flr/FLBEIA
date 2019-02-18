#----------------------------------------------------------------------------
#' Function to generate an FLStock from an FLBiol and FLFleet 
#'  
#'
#' @param biol An FLBiolobject
#' @param fleets An FLFleetExt object
#' 
#' @return An FLStock object with abundance and biological information from biol argument 
#'         and exploitation levels from fleets argument
#' 
#' 
#' @examples
#'\dontrun{
#'
#' data(res_flbeia)
#' 
#' biol   <- oneRes$biols$stk1
#' fleets <- oneRes$fleets
#' 
#' stk <- biolfleets2flstock( biol=biol, fleets=fleets)
#' 
#'}
#
# Sonia Sanchez
# 2019/02/04
#----------------------------------------------------------------------------

biolfleets2flstock <- function( biol, fleets){ 
  
  # checkings
  
  if( class(biol)!="FLBiol") 
    stop("biol object must be of class 'FLBiol'")
  
  if( class(fleets)!="FLFleetsExt") 
    stop("fleet object must be of class 'FLFleetsExt'")
  
  
  # Perfect observation
  res <- perfectObs( biol, fleets, covars=NULL, obs.ctrl=NULL, 
                     year=dim(biol@n)[2])
  
  return(res)
  
}

# Function to generate an FLStock from an FLBiol and FLFleet 
# (defined function is only valid for simple cases with only one unique stock with one unit, season, area and iter)
#
# biolfleets2flstock <- function( biol, fleets){ 
#   # checkings
#   if( class(biol)!="FLBiol") 
#     stop("biol object must be of class 'FLBiol'")
#   if( class(fleets)!="FLFleetsExt") 
#     stop("fleet object must be of class 'FLFleetsExt'")
#   if (sum(dim(biol@n)[3:6]>1)) 
#     stop("unit, season, area or iter has dimension higher than 1")
#   # Dimensions
#   na <- dim(biol@n)[1]
#   yr <- dim(biol@n)[2]
#   # Creating and filling FLStock
#   res <- as(biol, "FLStock")[,1:(yr-1),1,1]                          # stock.n, stock.wt, m, mat
#   stock(res) <- quantSums(res@stock.n*res@stock.wt)                  # stock
#   landings.n(res) <- apply(landStock(fleets, name(biol)), c(1:2,6),sum)[,1:(yr-1),]  # landings.n   
#   discards.n(res) <- apply(discStock(fleets, name(biol)), c(1:2,6),sum)[,1:(yr-1),]  # discards.n
#   catch.n(res)    <- res@discards.n + res@landings.n                 # catch.n
#   res@landings.wt <- res@discards.wt <- res@catch.wt <- res@stock.wt # weights at age
#   landings(res)   <- quantSums(res@landings.n*res@landings.wt)       # landings
#   discards(res)   <- quantSums(res@discards.n*res@discards.wt)       # discards
#   catch(res)      <- res@landings + res@discards                     # catch
#   # harvest
#   if(na == 1){ # * if age structured calculate it from 'n'.
#     harvest(res) <- (res@stock.n*res@stock.wt) * (1/res@catch)
#     units(res@harvest) <- "hr"
#   } else {      # * if biomass dyn => assume C = q*E*B => C = F*B and F = C/B.
#     units(res@harvest) <- "f"
#     res@harvest[-c((na-1):na),] <- log(biol@n[-c((na-1):na),-yr]/biol@n[-c(1,na),-1]) - 
#       res@m[-c((na-1):na),]
#     res@harvest[-c((na-1),na),] <- ifelse( is.na(res@harvest[-c((na-1),na),]) | 
#                                              as.numeric(res@harvest[-c((na-1),na),])<0, 
#                                            0, res@harvest[-c((na-1),na),])
#     # for older ages: optimisation
#     n. <- array(res@stock.n[drop=T], dim = c(na,yr-1)) # [na,ny]
#     m. <- array(res@m[drop=T], dim = c(na,yr-1))       # [na,ny]
#     c. <- array(res@catch.n[drop=T], dim = c(na,yr-1)) # [na,ny]
#     
#     fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
#     
#     for(y in 1:(yr-1)) for(a in (na-1):na) {
#       if(is.na(n.[a,y])) { res@harvest[a,y,,,,]<-0
#       } else {
#         if(n.[a,y] < c.[a,y]) { res@harvest[a,y,,,,] <- 10
#         } else {
#           zz <- try( ifelse(n.[a,y] == 0 | c.[a,y] == 0, 0,
#                             uniroot(fobj, lower = 1e-300, upper = 1e6, n = n.[a,y], m=m.[a,y], c = c.[a,y])$root), 
#                      silent = TRUE)
#           res@harvest[a,y,,,,] <- ifelse(is.numeric(zz), zz, res@harvest[(na-1)-1,y,,,,] )
#         }
#       }
#     }
#     harvest.spwn(res) <- m.spwn(res)
#   }
#   return(res)
# }


