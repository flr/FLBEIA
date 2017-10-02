# FLBEIA
- Version: 1.15.2
- Date: 2017-10-02
- Author: Dorleta GARCIA <dgarcia@azti.es>
- Maintainer: Dorleta GARCIA, AZTI.
- Repository: <https://github.com/flr/FLBEIA/>
- Bug reports: <https://github.com/flr/FLBEIA/issues>

## Overview
Provides a simulation toolbox which facilitates the development of bio-economic impact assessments of fisheries management strategies. It is built under a Management Strategy Evaluation framework using FLR. The simulation is divided in two *worlds*: the operating model (OM, i.e. the real world) and the management procedure model (MP, i.e. the perceived world). The OM is itself divided in 3 components: biological, fleets and environmental and economic covariates. The MP is also divided in 3 components: data, the perceived system and advice.

To install this package, start R and enter:

	install.packages(c("plyr", "ggplot2", "nloptr", "mvtnorm", "triangle"))
	install.packages("FLBEIA", repos="http://flr-project.org/R")

## Documentation
- [Help pages](http://flr-project.org/FLBEIA)
- Vignette

## Citation

To cite this package, start R and enter

	library(FLBEIA)
	citation("FLBEIA")

## License
Copyright (c) 2010-2017 AZTI. Released under the [GPL v2](http://www.gnu.org/licenses/gpl-2.0.html).

## Contact
You are welcome to:

- Submit suggestions and bug-reports at: <https://github.com/flr/FLBEIA/issues>
- Send a pull request on: <https://github.com/flr/FLBEIA/>
- Compose a friendly e-mail to: <dgarcia AT azti.es> or <flbeia AT azti.es>
