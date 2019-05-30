# The R Package: DendroSync #

A Set of Tools for Calculating Spatial Synchrony Between Tree-Ring Chronologies.

### General Description ###

* Dendrosync is an R library that provides functions for the calculation and plotting of synchrony in the tree growth from tree-ring width chronologies (TRW index). It combines variance-covariance (VCOV) mixed modelling with functions that quantify the degree to which TRW chronologies contain a common temporal signal. It also implements temporal trends in spatial synchrony using a moving window. These methods can also be used with other kind of eacological variables that have temporal autocorrelation corrected. 

* Version 0.1.3

* Depends: R (>= 3.1.2), nlme, ggplot2

### Package Contents ###

* The package contains different functions to evaluate synchrony in tree growth from tree-ring width chronologies. Each function is described in help files with clear examples.
For a detailed description of the package visit:
[josugalday/r-stuff](https://sites.google.com/view/josugalday/r-stuff)

* Main functions:

1. `dendro.varcov`: it calculates variance-covariance (VCOV) mixed models

2. `mod.table`: it provides a table to compare fitted variance-covariance (VCOV) mixed models by AIC, AICc, BIC

3. `sync`: it calculates spatial synchrony from fitted variance-covariance mixed models  

4. `sync.plot`: it creates dot plots of within- and between-group synchrony

5. `sync.trend`: it calculates temporal trends of spatial sinchrony

6. `sync.trend.plot`: it creates a line chart showing temporal trends of spatial synchrony 

### Contact ###

* For any problem please contact with Josu G Alday (josucham@gmail.com) or Victor Resco (v.rescodedios@gmail.com)

* Citation: 

Alday JG; Shestakova TA, Resco de Dios V, Voltas J. (2018) DendroSync: An R package to unravel synchrony patterns in tree-ring networks. Dendrochronologia 47: 17-22.