# WhyWhere

This package is an R implemention of the WhyWhere data mining algorithm. Developed for biodiversity modelling with big data issues in mind, the approach is simple and rigorous, efficient to compute, and provides accuracy equal to the best alternative approaches.  It represents the way forward in utilize large spatial environmental data sets in a robust predictive and explanatory capacity.  

Current release: v0.0.2 on GitHub.

Introduction, installation, documentation, benchmarks etc to be found: [HOMEPAGE](https://github.com/davids99us/whywhere/wiki)

To install:

    >install.packages("devtools")
    >library(devtools)
    >install_github("davids99us/whywhere")
    
Global data file called global.tar.kz can be downloaded by ftp from landshape.org contains around 1000 global spatial data layers in a wide range of themes and resolutions.  This data is 179MBytes in (compressed) size.
    
A vignette is available in the vignette directory.

See an example that tests the package by running test().


##NEW FEATURES

    Now passes R CMD check but still some warnings.
    
    