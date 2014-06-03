library(KEGGBioCycAPI)
library(foreach)
library(doMC)
registerDoMC(4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create the whole a KEGG and BioCyc List.
# Afer that, two dataset 'biocycSpe.RData' and 'bioKEGGSpe.RData' will be created.
# ATTENTION: The KEGG and BioCyc database may update frequently, so please re-run the code below to retrieve the updated species information. Last update: July 06, 2014.

# whole KEGG species list
wKEGGSpe <- getSpePhylo(whole = TRUE)
wNCBISpe <- KEGG2Tax(wKEGGSpe[, 2])
wNCBISpe <- wNCBISpe[order(names(wNCBISpe))]
wNCBISpe <- wNCBISpe[rank(wKEGGSpe[,2])]
wKEGGSpe <- cbind(wKEGGSpe, wNCBISpe)
colnames(wKEGGSpe)[5] <- 'TaxonomyID'
save(wKEGGSpe, file = 'bioKEGGSpe.RData')

# whole BioCyc species list
wBioCyc <- getPhyloCyc(whole = TRUE)
wNCBISpe <- cyc2Tax(wBioCyc[, 1])
wNCBISpe <- wNCBISpe[order(names(wNCBISpe))]
wNCBISpe <- wNCBISpe[rank(wBioCyc[,1])]
wBioCyc <- cbind(wBioCyc, wNCBISpe)
colnames(wKEGGSpe)[4] <- 'TaxonomyID'
save(wBioCyc, file = 'biocycSpe.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load whole KEGG and BioCyc species information
load('biocycSpe.RData')
load('bioKEGGSpe.RData')
