library(KEGGBioCycAPI)
library(foreach)
library(doMC)
library(XML)
registerDoMC(4)


# load whole KEGG and BioCyc species information
load('biocycSpe.RData')
load('bioKEGGSpe.RData')
