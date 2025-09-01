
# isoTB: Stable Isotope Geochemistry Tool Box

## Overview

The isoTB package provides multiple tools to help process empirical or
model stable isotope geochemistry data for environmental science
studies. The isoTB package provides simple conversion functions to
convert between different isotopic composition notations (*e.g.* isotope
ratio, delta notation, fractional abundance), isotope mixing models,
isotope mixing model solvers and functions to calculate the isotope
fractionation of sources.

The isoTB package also provides very useful functions for enriched
isotope tracer studies enabling (i) to calculate the error induced by
neglecting isotope fractionation processes in tracing studies and (ii)
to calculate the minimum isotope enrichment of the tracer source in a
system to ensure that the aforementioned error is below a user-defined
threshold.

isoTB can be used to model complex systems for both natural variations
of isotope ratios and enriched isotope tracer studies by creating the
different reservoirs in the system and defining the interactions
(isotope mixing, isotope fractionation) between the different
reservoirs.

## Installation

You can install the development version of isoTB from GitHub like so:

``` r
# Using the pak package to install directly from GitHub
# install.packages("pak")
pak::pak("gregoryVDH/isoTB") 

# Using the devtools package to install directly from GitHub
# install.packages("devtools")
devtools::install_github("gregoryVDH/isoTB") 
```
