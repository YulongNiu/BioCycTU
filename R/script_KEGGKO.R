library(KEGGBioCycAPI)
library(foreach)
library(doMC)
registerDoMC(4)

#################### Get whole KEGG and BioCyc species lists #################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create the whole a KEGG and BioCyc List.
# Afer that, two dataset 'wBiocycSpe.RData' and 'wKEGGSpe.RData' will be created.
# ATTENTION: The KEGG and BioCyc database may update frequently, so please re-run the code below to retrieve the updated species information. Last update: July 06, 2014.

# whole KEGG species list
wKEGGSpe <- getSpePhylo(whole = TRUE)
wNCBISpe <- KEGG2Tax(wKEGGSpe[, 2])
wNCBISpe <- wNCBISpe[order(names(wNCBISpe))]
wNCBISpe <- wNCBISpe[rank(wKEGGSpe[,2])]
wKEGGSpe <- cbind(wKEGGSpe, wNCBISpe)
colnames(wKEGGSpe)[5] <- 'TaxonomyID'
save(wKEGGSpe, file = 'wKEGGSpe.RData')

# whole BioCyc species list
wBiocycSpe <- getPhyloCyc(whole = TRUE)
wNCBISpe <- cyc2Tax(wBiocycSpe[, 1])
wNCBISpe <- wNCBISpe[order(names(wNCBISpe))]
wNCBISpe <- wNCBISpe[rank(wBiocycSpe[,1])]
wBiocycSpe <- cbind(wBiocycSpe, wNCBISpe)
colnames(wBiocycSpe)[4] <- 'TaxonomyID'
save(wBiocycSpe, file = 'wBiocycSpe.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

########### Preprocess the species list from KEGG and BioCyc ############

# load whole KEGG and BioCyc species information
load('wKEGGSpe.RData')
load('wBiocycSpe.RData')

# preprocess the BioCyc database. Remove the 'TaxonomyID' with ""
wBiocycSpe <- wBiocycSpe[!(wBiocycSpe[, 'TaxonomyID'] %in% ''), ]
# rm the BioCycID with the suffix with 'HMP'
## taxDup <- wBiocycSpe[duplicated(wBiocycSpe[, 'TaxonomyID']), 'TaxonomyID']
## test1 <- wBiocycSpe[wBiocycSpe[, 'TaxonomyID'] %in% taxDup, ]
## test1 <- test1[order(test1[, 'TaxonomyID']), ]
## test2 <- wKEGGSpe[wKEGGSpe[, 'TaxonomyID'] %in% taxDup, ]
## test2 <- test2[order(test2[, 'TaxonomyID']), ]
## write.csv(test1, 'test1.csv')
## write.csv(test2, 'test2.csv')
reSpe <- c('SSPU546271-HMP', 'PDEN908937-HMP', 'LPAR537973-HMP', 'FSP469604-HMP', '', 'ECOL566546-HMP', 'EFAE333849-HMP', 'CEFF196164-HMP', 'CAUR548476-HMP', 'BABO359391')
wBiocycSpe <- wBiocycSpe[!(wBiocycSpe[, 'BioCycID'] %in% reSpe), ]



# preprocess the KEGG database. Remove the duplicated 'TaxonomyID'
## taxDup <- wKEGGSpe[duplicated(wKEGGSpe[, 'TaxonomyID']), 'TaxonomyID']
## tmp1 <- wKEGGSpe[wKEGGSpe[, 'TaxonomyID'] %in% taxDup, ]
## tmp1 <- tmp1[order(tmp1[, 'TaxonomyID']), ]
## tmp2 <- wBiocycSpe[wBiocycSpe[, 'TaxonomyID'] %in% taxDup, ]
## tmp2 <- tmp2[order(tmp2[, 'TaxonomyID']), ]
## test1 <- foreach(i = 1:nrow(tmp1)) %dopar% {
##   t1 <- getKEGGPathGenes(tmp1[i, 2])
##   t2 <- t1[names(t1) %in% paste('path:', tmp1[i, 2], '00190', sep = '')]
##   t2 <- t2[[1]]
## }
## names(test1) <- tmp1[, 2]
## write.csv(tmp1, 'tmp1.csv')
## write.csv(tmp2, 'tmp2.csv')
reSpe <- c('scy', 'css', 'syz', 'syy', 'ply', 'tae', 'plr', 'cgb', 'bms', 'tmm', 'tmi', 'tpw', 'msg', 'eru', 'gdj', 'lpu', 'bld', 'chb', 'vca', 'vcr', 'bafz', 'blon', 'bmu', 'osa', 'dosa', 'mtur', 'ebe', 'mre', 'oca', 'edj', 'ell', 'lrh', 'cty', 'fsc', 'ekf', 'hte', 'rai', 'amn', 'mtv', 'heo', 'lpo')
wKEGGSpe <- wKEGGSpe[!(wKEGGSpe[, 'KEGGID'] %in% reSpe), ]

#################### transfer KEGG gene ID to BioCyc gene ID #################
# get A1Ao ATP synthase KO
AATPKO <- c('K02117', 'K02118', 'K02119', 'K02120', 'K02121', 'K02122', 'K02107', 'K02123', 'K02124')
AATPKONames <- c('A', 'B', 'C', 'D', 'E', 'F', 'H', 'I', 'K')
AATPKObySpe <- KEGGSpeKO(AATPKO, AATPKONames)

# merge species lists
mergeList <- merge(wBiocycSpe, wKEGGSpe, by.x = 'TaxonomyID', by.y = 'TaxonomyID')

# transfer to BioCyc species ID
commSpeIDs <- intersect(names(AATPKObySpe), mergeList[, 'KEGGID'])
commAATPKO <- AATPKObySpe[names(AATPKObySpe) %in% commSpeIDs]
commList <- mergeList[mergeList[, 'KEGGID'] %in% commSpeIDs, ]
commList <- commList[order(commList[, 'KEGGID']), ]
commList <- commList[rank(names(commAATPKO)), ]
names(commAATPKO) <- commList[, 'BioCycID']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# transfer KEGG gene IDs to BioCyc gene IDs
# BioCyc geneID. Cut the whole length with internal 4
cutMat <- CutSeqEqu(length(commAATPKO), 4)

for (j in 59:ncol(cutMat)) {
  ATPKOCycPart <- foreach (i = cutMat[1, j]:cutMat[2, j]) %dopar% {
    # may have NA
    print(paste('It is running ', j, ' in cutMat. ', 'It is running ', i, ' with the name of ', names(commAATPKO)[i], '.', sep = ''))
    iniVal <- commAATPKO[[i]]
    iniVal[!is.na(iniVal)] <- sapply(iniVal[!is.na(iniVal)], KEGGID2CycID, speKEGGID = commList[i, 'KEGGID'], speCycID = names(commAATPKO)[i])
    return(iniVal)
  }

  names(ATPKOCycPart) <- names(commAATPKO)[cutMat[1, j]:cutMat[2, j]]
  save(ATPKOCycPart, file = paste('ATPKOCycPart', cutMat[1, j], '_', cutMat[2, j], '.RData', sep = ''))
  Sys.sleep(30)
}

# read all the BioCyc genes
partFile <- dir(pattern="^ATPKOCycPart.*")
ATPCycList <- foreach(i = 1:length(partFile), .combine = append) %dopar% {
  load(partFile[i])
  return(ATPKOCycPart)
}

# deal with the list
## sink(file = 'tmp1.txt')
## ATPCycList[sapply(ATPCycList, is.list)]
## sink()
ATPCycList$CSTI499177$A <- 'GJE9-2386'
ATPCycList$CSTI499177$B <- 'GJE9-2385'
ATPCycList$HELO768066$A <- 'GJEE-1718'
ATPCycList$HELO768066$B <- 'GJEE-1724'
ATPCycList$HELO768066$D <- 'GJEE-1723'
ATPCycList$HELO768066$I <- 'GJEE-1722'
ATPCycList$LLON661367$A <- 'GJAR-338'
ATPCycList$LLON661367$B <- 'GJAR-332'
ATPCycList$MACE188937$C <- 'GI2O-4209'
ATPCycList$MACE188937$K <- 'GI2O-4207'
ATPCycList$`MALC1091494-WGS`$A <- 'GSQK-1672'
ATPCycList$`MALC1091494-WGS`$B <- 'GSQK-1666'
ATPCycList$PAMO264201$B <- 'GH0M-1715'
ATPCycList$BHYO565034$A <- 'GJI7-1902'
ATPCycList$BHYO565034$B <- 'GJI7-1903'
ATPCycList$BHYO565034$E <- 'GJI7-1900'
ATPCycList$BPIL1133568$B <- 'GL9N-312'
ATPCycList$BPIL759914$B <- 'GHZ5-1034'
ATPCycList$CNOV386415$B <- 'GH98-780'
ATPCycList <- sapply(ATPCycList, unlist)

# deal with "0"
## test1 <- sapply(ATPCycList, function(x) {
##   if(sum(x == "0", na.rm = TRUE) > 0) {
##     return(TRUE)
##   } else {return(FALSE)}
## })
## sink(file = 'tmp1.txt')
## ATPCycList[test1]
## sink()
ATPCycList$`FACI333146-WGS`['A'] <- 'GSNJ-959'
ATPCycList$`FACI333146-WGS`['B'] <- 'GSNJ-958'
ATPCycList$`FACI333146-WGS`['C'] <- 'GSNJ-961'
ATPCycList$`FACI333146-WGS`['D'] <- 'GSNJ-957'
ATPCycList$`FACI333146-WGS`['E'] <- 'GSNJ-962'
ATPCycList$`FACI333146-WGS`['F'] <- 'GSNJ-960'
ATPCycList$`FACI333146-WGS`['H'] <- 'GSNJ-955'
ATPCycList$`FACI333146-WGS`['I'] <- 'GSNJ-954'
ATPCycList$`FACI333146-WGS`['K'] <- 'GSNJ-963'
ATPCycList$AXYL698758['C1'] <- 'GL81-572'
ATPCycList$AXYL698758['F1'] <- 'GL81-575'
ATPCycList$AXYL698758['I1'] <- 'GL81-573'
ATPCycList$AXYL698758['K1'] <- 'GL81-574'
ATPCycList$`TNIT1255043-WGS`['B'] <- 'GSYW-3805'
ATPCycList$`TNIT1255043-WGS`['D'] <- 'GSYW-3804'
ATPCycList$`TNIT1255043-WGS`['E'] <- 'GSYW-3807'
ATPCycList$`TNIT1255043-WGS`['I'] <- 'GSYW-3810'
ATPCycList$CDIF272563['A'] <- 'GJFE-3189'
ATPCycList$CDIF272563['B'] <- 'GJFE-3188'
ATPCycList$CDIF272563['C'] <- 'GJFE-3191'
ATPCycList$CDIF272563['D'] <- 'GJFE-3187'
ATPCycList$CDIF272563['E'] <- 'GJFE-3192'
ATPCycList$CDIF272563['F'] <- 'GJFE-3190'
ATPCycList$CDIF272563['H'] <- 'GJFE-3195'
ATPCycList$CDIF272563['I'] <- 'GJFE-3194'
ATPCycList$CDIF272563['K'] <- 'GJFE-3193'
ATPCycList$CMUR243161['E'] <- 'GHYU-604'
ATPCycList$CMUR243161['I'] <- 'GHYU-599'
ATPCycList$MBOU1201294['A'] <- 'GLGD-173'
ATPCycList$MBOU1201294['F'] <- 'GLGD-171'
ATPCycList$MBOU1201294['K'] <- 'GLGD-168'
# delete TOSH751945
ATPCycList <- ATPCycList[!(names(ATPCycList) %in% 'TOSH751945')]
save(ATPCycList, file = 'ATPCycList.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
