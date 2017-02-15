
# FLBEIA
- Version: 1.15
- Date: 2016-02-16
- Author: Dorleta GARCIA <dgarcia AT azti.es>
- Maintainer: Dorleta GARCIA <dgarcia AT azti.es>
- Repository: <https://github.com/flr/FLBEIA/>
- Bug reports: <https://github.com/flr/FLBEIA/issues>

## Overview
`FLBEIA` provides a simulation toolbox which facilitates the development of bio-economic impact assessments of fisheries management strategies. The objective of the model is to facilitate the Bio-Economic Impact Assessment of a great Variety of Management strategies. The model is multistock, multifleet and seasonal. The simulation is divided in 2 main parts, the Operating Model (OM) and the Management Procedure Model (MPM). At the same time the OM is divided in 3 components, the biological component, which simulates the stocks dynamics, the fleets component, that simulates the fleets' dynamics and the covariables component, that simulates the dynamic of any variable that is not included in biological and fleets component. The MPM is also divided in 3 components, the observation model, which simulates the data collection, the assessment model, which uses the simulated data to produce a 'observed population' and the advice model which produces the advice based on the output of the assessment model.

To install this package, start R and enter:
  
  install.packages(c("plyr", "ggplot2", "nloptr", "mvtnorm", "triangle"))
	install.packages("FLBEIA", repos="http://flr-project.org/R")

or download from the [FLBEIA releases page](https://github.com/flr/FLBEIA/releases/latest)

## Documentation
- [Help pages](http://flr-project.org/FLBEIA/Reference)

## Build Status
[![Travis Build Status](https://travis-ci.org/flr/FLBEIA.svg?branch=master)](https://travis-ci.org/flr/FLBEIA)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/flr/FLBEIA?branch=master&svg=true)](https://ci.appveyor.com/project/flr/FLBEIA)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/FLBEIA)](https://cran.r-project.org/package=FLBEIA)

## Releases
- [All release](https://github.com/flr/FLBEIA/releases/)

## License
Copyright (c) 2004-2015 The FLR Team. Released under the [GPL v2](http://www.gnu.org/licenses/gpl-2.0.html).

## Contact
You are welcome to:

- Submit suggestions and bug-reports at: <https://github.com/flr/FLBEIA/issues>
- Send a pull request on: <https://github.com/flr/FLBEIA/>
- Compose a friendly e-mail to: <dgarcia AT azti.es>
