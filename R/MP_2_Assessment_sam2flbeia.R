###########################################################
#
#               SAM 
#
# Developed by Ibrahim Umar (IMR Norway) 20/11/2018
# NOTE: Install FLSAM (v2):
# > devtools::install_github("flr/FLSAM")
#..........................................................


sam2flbeia <- function(stock, indices, control=control, covars=covars){

  stock@landings.n[stock@landings.n == 0] <- 0.1
  stock@catch.n <- stock@landings.n + stock@discards.n

  # We have to specify the type of the indices
  for(i in 1:length(indices))
	indices[[i]]@type <- control$indices.type

  stock.ctrl <- FLSAM.control(stock, indices)
  fit <- FLSAM(stock, indices, stock.ctrl, return.fit = TRUE)
  stock.sam <- SAM2FLR(fit, stock.ctrl)
  stock <- stock + stock.sam

  return(list(stock=stock,covars=covars))

}



