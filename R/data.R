#-------------------------------------------------------------------------------  
# data description: 
# description of available datasets in FLBEIA library
# 
# Created: Sonia Sanchez -  2018-07-05
# Changed: 
#------------------------------------------------------------------------------- 

# data.R - description of datasets
# FLBEIA/R/data.r

# Copyright: AZTI, 2018
# Author: FLBEIA team (AZTI) (<flbeia@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

#' FLBEIA datasets
#'
#' Example datasets for the classes defined in FLBEIA.
#'
#' \itemize{
#' 
#'    \item{\code{one}}{: A dataset for running FLBEIA.
#'                         Example with one stock (age-structured), one fleet, annual steps (one season) and one iteration.}
#'     \itemize{
#'       \item{\code{oneBio} (\code{\link{FLBiols}})}
#'            {: Biological information on the stock. In this case, the stock is age-structured.}
#'       \item{\code{oneSR} (\code{\link{FLSRsim}})}
#'            {: Stock-recruitment model for the stock in \code{oneBio} object.}
#       \item{\code{BDs} (\code{\link{FLSRsim}})}
#            {: Biomass dynamic model for the stock in biomass in \code{oneBio} object.}
#'       \item{\code{oneFl} (\code{\link{FLFleetsExt}})}
#'            {: Information on the fleet and the metier considered.}
#'       \item{\code{oneCv} (\code{list} of \code{\link{FLQuants}})}
#'            {: Covariates information. In this case all are economic indicators.}
#'       \item{\code{oneIndAge} and \code{oneIndBio} (\code{list} of \code{\link{FLIndices}})}
#'            {: Indices, if avalable, for the different stocks in \code{oneBio}.
#'               Where \code{oneIndAge} and \code{oneIndBio} are indices with estimates in numbers at age and total biomass, respectively.}
#'       \item{\code{oneAdv} (\code{list})}
#'            {: Information on TAC and quota share.}
#'       \item{\code{oneMainC} (\code{list})}
#'            {: Settings to control the main function \code{\link{FLBEIA}}. The simulation years.}
#'       \item{\code{oneBioC} (\code{list})}
#'            {: Settings to control the biological operating model for the stock in \code{oneBio} object.}
#'       \item{\code{oneFlC} (\code{list})}
#'            {: Settings to control the fleet operating model for the fleet in \code{oneFl} object.}
#'       \item{\code{oneCvC} (\code{list})}
#'            {: Settings to control the covar operating model for each covariate in \code{covars} object.
#'               In this case all the covariates are fixed.}
#'       \item{\code{oneObsC}, \code{oneObsCIndAge} and \code{oneObsCIndBio} (\code{list})}
#'            {: Settings to control the observation model for the stock in \code{oneBio} object.
#'               In \code{oneObsC} the stock is observed without error and there are no indices available.
#'               Alternative control settings are available for cases when indices are observed, \code{oneObsCIndAge} and \code{oneObsCIndBio}.}
#'       \item{\code{oneAssC} (\code{list})}
#'            {: Settings to control the assessment model for each stock in \code{oneBio} object.
#'               In this case, no assessment is carried out.}
#'       \item{\code{oneAdvC} (\code{list})}
#'            {: Settings to control the advice model for the stock in \code{oneBio} object.}
#'     }
#'     
#'    \item{\code{oneIt}}{: A dataset for running FLBEIA.
#'                          Same as \code{one} dataset, but with three iterations.}
#'     
#'    \item{\code{multi}}{: A dataset for running FLBEIA.
#'                          Example with two stocks (one age-structured and the other in biomass), 
#'                          two fleets (with 2 metiers each), four seasons and one iteration.}
#'     \itemize{
#'       \item{\code{multiBio} (\code{\link{FLBiols}})}
#'            {: Biological information on the stocks (one is age-structured and the other one in total biomass).}
#'       \item{\code{multiSR} (\code{\link{FLSRsim}})}
#'            {: Stock-recruitment models for the age-structured stock in \code{multiBio} object.
#'               In this case a Beverton-Holt (\code{\link{bevholt}}) is selected.}
#'       \item{\code{multiBD} (\code{\link{FLSRsim}})}
#'            {: Biomass dynamic model for the stock in biomass in \code{multiBio} object.
#'               In this case Pella-Tomlinson model (\code{PellaTom}) is selected.}
#'       \item{\code{multiFlC} (\code{\link{FLFleetsExt}})}
#'            {: Information on the fleets and metiers. 
#'               In this case there are two fleets, each one with two metiers, all of them capturing both stocks in \code{multiBio} object.}
#'       \item{\code{multiCv} (\code{list} of \code{\link{FLQuants}})}
#'            {: Covariates information. In this case all are economic indicators.}
#       \item{\code{indices} (\code{list} of \code{\link{FLIndices}})}
#            {: Indices, if avalable, for the different stocks in \code{multiBio}.}
#'       \item{\code{multiAdv} (\code{list})}
#'            {: Information on TAC and quota share.}
#'       \item{\code{multiMainC} (\code{list})}
#'            {: Settings to control the main function \code{\link{FLBEIA}}.}
#'       \item{\code{multiBioC} (\code{list})}
#'            {: Settings to control the biological operating model for each stock in \code{multiBio} object.}
#'       \item{\code{multiFlC} (\code{list})}
#'            {: Settings to control the fleet operating model for each fleet in \code{fleets} object.}
#'       \item{\code{multiCvC} (\code{list})}
#'            {: Settings to control the covar operating model for each covariate in \code{covars} object. 
#'               In this case all the covariates are fixed.}
#'       \item{\code{multiObsC} (\code{list})}
#'            {: Settings to control the observation model for each stock in \code{multiBio} object. 
#'               In this case the stock is observed without error and there are no indices available. }
#'       \item{\code{multiAssC} (\code{list})}
#'            {: Settings to control the assessment model for each stock in \code{multiBio} object. 
#'               In this case, no assessment is carried out.}
#'       \item{\code{multiAdvC} (\code{list})}
#'            {: Settings to control the advice model for each stock in \code{multiBio} object.}
#'     }
#'     
#'    \item{\code{res_flbeia}}{: A dataset with the outputs of FLBEIA runs given as input the different datasets (one, oneIt and multi).}
#'     \itemize{
#'       \item{\code{oneRes} (\code{list})}
#'            {: Output of the FLBEIA function, given as input the data in the \code{one} dataset.}
#'       \item{\code{onIteRes} (\code{list})}
#'            {: Output of the FLBEIA function, given as input the data in the \code{oneIt} dataset.}
#'       \item{\code{multiRes} (\code{list})}
#'            {: Output of the FLBEIA function, given as input the data in the \code{multi} dataset.}
#'     }
#'     
#'    \item{\code{mur}}{: A dataset for Stripped Red Mullet in the Bay of Biscay.
#'                        Information on catch and abundance indices from Evohe survey.}
#'     \itemize{
#'       \item{\code{catch} (\code{\link{data.frame}})}
#'            {: The total catch time series data by area from WGBIE report (ICES, 2017).
#'               Total catch data is available since 1975. In 1999 France did not report any data. 
#'               As France is the main contributor to the total catch, the 1999 catch data was not included in the analysis.}
#'       \item{\code{evhoe} (\code{\link{data.frame}})}
#'            {: EVHOE abundance index time series, provided by Ifremer.
#'               The abundance index is available since 1997 and provides an estimation of the biomass together with a coefficient of variation.}
#'     }
#' }
#'
#' Datasets can be loaded by issuing the \code{data} command, like in:
#' \code{data(one)}.
#' 
#' All available datasets can be checked by: \code{data(package='FLBEIA')}.
#'
#' @name datasets
#' @aliases one oneIt multi res_flbeia mur
#' @seealso \link{FLBEIA}, \linkS4class{FLBiols}, \linkS4class{FLFleetsExt}, \linkS4class{FLSRSim},
#' \linkS4class{FLIndices}, \linkS4class{FLQuant}
#' 
# @references
#' @keywords datasets
#' @references 
#'   ICES, 2017
#' 
#' @examples
#'
#' data(one)
#' data(res_flbeia)
#'

NULL
