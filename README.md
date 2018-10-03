# Overview

**grase** is a package that estimates the sequencing effort required to achieve a goal of MAG construction. This is based on a model presented in [Royalty and Steen \(2018\)](https://doi.org/10.1101/356840).

Note that this method simply estimates coverage of genomes in a mixed community; it doesn't take into account efficiency of assembly or accuracy of binning algorithms. It therefore places a lower bound on the sequencing effort that is required to avhieve a given coverage.

# Use

**grase** can be run over the web or locally. Running over the web is the easiest: navigate to [adsteen.shinyapps.io/grase](http://adsteen.shinyapps.io/grase). Alternately, you can run on your own computer using the following script:

```
library(shiny)
runGitHub("adsteen/grase/grase.gz)
```

**grase** relies on the following packages, which must be installed for it to run locally: **mgcv**, **shiny**, **shinyWidgets**, **dplyr**, **ggplot2**, **scales**.**
