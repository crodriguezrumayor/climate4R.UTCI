# climate4R.UTCI

## What is `climate4R.UTCI`?

A R package for computing the Universal Thermal Climate Index (UTCI) directly from climate data within the `climate4R` framework.

`climate4R.UTCI` is a wrapper of the [`UTCIr`](https://github.com/ECA-D/UTCIr) R package by ECA&D, which implements the UTCI using a vectorised Fortran backend. This wrapper adapts that function for seamless integration with the **climate4R** data structures, providing support for parallel computing.

## Installation

The recommended procedure for installing the package is using the `remotes` package:

```R
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("crodriguezrumayor/climate4R.UTCI")
```

Note that the following dependencies need to be installed beforehand:

```R
remotes::install_github("SantanderMetGroup/transformeR")
remotes::install_github("SantanderMetGroup/convertR")
remotes::install_github("ECA-D/UTCIr")
```

## Usage

```R
library(climate4R.UTCI)

data("ERA5_day_t2mx", package = "climate4R.UTCI")
data("ERA5_day_t2mn", package = "climate4R.UTCI")
data("ERA5_day_hurs", package = "climate4R.UTCI")
data("ERA5_day_sfcwind", package = "climate4R.UTCI")
data("ERA5_day_ssrd", package = "climate4R.UTCI")
data("ERA5_lsm", package = "climate4R.UTCI")       # optional land-sea mask

utci <- utciGrid(tmax = ERA5_day_t2mx, tmin = ERA5_day_t2mn,
                 hurs = ERA5_day_hurs, wind = ERA5_day_sfcwind,
                 radiation = ERA5_day_ssrd, mask = ERA5_lsm)
```

## References

Bröde, P. et al. (2012). Deriving the operational procedure for the Universal Thermal Climate Index (UTCI). International Journal of Biometeorology, 56, 481-494. https://doi.org/10.1007/s00484-011-0454-1
