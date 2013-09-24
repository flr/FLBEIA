# FLBEIA

**Bio-Economic Impact Assessment of Management strategies using FLR**

Version: 1.0

`FLBEIA` provides a simulation toolbox which facilitates the development of bio-economic impact assessments of fisheries management strategies. It is built under a Management Strategy Evaluation framework using FLR. The simulation is divided in two *worlds*: the operating model (OM, i.e. the real world) and the management procedure model (MP, i.e. the perceived world). The OM is itself divided in 3 components: biological, fleets and environmental and economic covariates. The MP is also divided in 3 components: data, the 
perceived system and advice.

Author: Dorleta GARCIA <dgarcia@azti.es>

Maintainer: Dorleta GARCIA, AZTI Tecnalia.


## Installation

**TO INSTALL** this package, start R and enter:

	install.packages(repos="http://flr-project.org/Rdevel")

or, to install the development version

	library(devtools)
	install_github("FLBEIA", "flr")

**TO CITE** this package, start R and enter

	library(FLBEIA)
	citation("FLBEIA")

## Documentation

## Details

- Version: 1.3.20130719
- Date: 2013-07-19
- License: GPL-2
- Depends: R(>= 2.15.0), FLCore
- Imports: FLCore
- URL: <https://github.com/flr/FLBEIA>

### Classes
  - FLBDsim
  - FLSRsim
  - FLCatchesExt
  - FLCatchExt
  - FLFleetExt
  - FLFleetsExt
  - FLMetierExt
  - FLMetiersExt


### Methods
  - addFLCatch
  - catchNames
  - catchNames<-
  - FLCatchesExt
  - FLCatchExt
  - FLFleetExt
  - FLFleetsExt
  - FLMetierExt
  - FLMetiersExt
  - is.FLBDSRsims
  - is.FLCatchesExt
  - is.FLFleetsExt
  - is.FLMetiersExt

## Downloads
- Packaged Source [FLBEIA_1.3.20130719.tar.gz](http://flr-project.org/Rdevel/src/contrib/FLBEIA_1.3.20130719.tar.gz)
- Windows Binary [FLBEIA_1.3.20130719.zip](http://flr-project.org/Rdevel/bin/windows/contrib/3.0/FLBEIA_1.3.20130719.zip)

## Feedback
Please report any bug or issue in the [FLBEIA issue tracker](https://github.com/flr/FLBEIA/issues)
