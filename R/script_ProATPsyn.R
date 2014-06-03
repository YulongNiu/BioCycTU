#! /usr/bin/Rscript --vanilla

## subunit KOID
## beta K02112
## alpha K02111
## gamma K02115
## delta K02113
## epsilon K02114
## c K02110
## a K02108
## b K02109

CutSeq <- function(cutSeq){
  # USE: cut a vector based on the 'cutSeq'. This function is used to cut "full" seq. It means the sum of 'cutSeq' should be equal to the length of vector, which is used to cut.
  # INPUT: 'cutSeq' a sequence used to cut the vector.
  # OUTPUT: the index matrix

  # remove 0, because we cannot cut a sequence by the internal of 0.
  cutSeq <- cutSeq[cutSeq != 0]
  vecCutseq <- length(cutSeq)

  if (vecCutseq == 1) {
    headCut <- 1
    endCut <- cutSeq
  } else {
    # loopCutSeq is the circle of vecCutseq
    loopCutSeq <- list()
    for(i in 1:vecCutseq) {
      loopCutSeq[[i]] <- cutSeq[1:i]
    }

    loopSumCutSeq <-  sapply(loopCutSeq, sum)

    # the head and tail sequence
    headCut <- c(1,loopSumCutSeq[1:(vecCutseq-1)]+1)
    endCut <- loopSumCutSeq
  }

  cutMat <- matrix(c(headCut, endCut), 2, byrow=TRUE)

  return(cutMat)

}

CutSeqEqu <- function(vecLen, equNum){
  # USE: to cut a vector with equal internal.
  # INPUT: 'vecLen' the length of vector used to cut. 'equNum' the equal internal
  # OUTPUT: the index matrix.

  if (equNum > vecLen){
    # the internal is bigger than the length of vecLen. So we use the full vecLen.
    cutMat <- matrix(c(1, vecLen))
  } else {
    timeNum <- vecLen %/% equNum
    remainer <- vecLen %% equNum
    cutSeq <- c(rep(equNum, timeNum), remainer)
    cutMat <- CutSeq(cutSeq)
  }

  return(cutMat)
}


# load KEGG and BioCyc database
load('biocycSpe.RData')
load('bioKEGGSpe.RData')
load('KEGGCycAPI.RData')
library(foreach)
library(doMC)
library(XML)
registerDoMC(4)

# the common species
commSpe <- merge(biocycSpe, bioKEGGSpe, by.x = 'TaxonomyID', by.y = 'TaxonomyID', sort = FALSE)
commProSpe <- commSpe[grepl('Prokaryotes', commSpe[, 8]), ]
# remove duplicate BioCyc speID according to KEGG speID

commKEGGVec <- commProSpe[duplicated(commProSpe[, 6]), 6]
commKEGGMat <- commProSpe[commProSpe[, 6] %in% commKEGGVec, ]
delRowName <- c('553', '839', '1057', '1056', '1066', '1346', '1267', '1481', '1480', '1482', '1582', '1648', '1645', '1647', '1723', '1681', '1666', '1902', '2107', '1918', '2299', '2473')
commProSpe <- commProSpe[!(rownames(commProSpe) %in% delRowName), ]

# some species names changed
commProSpe[commProSpe[, 2] %in% 'BANT198094-WGS', 2] <- 'ANTHRA'

## # **select species used in phylogenetic profiling**
## phyloProSpe <- read.csv('wholeListFile.csv', row.names = 1)
## commProSpe <- commProSpe[commProSpe[, 6] %in% phyloProSpe[, 2], ]

# get KO list
KOVec <- c('K02111', 'K02112', 'K02115', 'K02113', 'K02114', 'K02110', 'K02108', 'K02109')
names(KOVec) <- c('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'c', 'a', 'b')
KOList <- vector('list', 8)
names(KOList) <- paste(names(KOVec), 'KO', sep = '')
for (i in 1:8) {
  KOMat <- getKEGGKO(KOVec[i])
  KOMat <- KOMat[KOMat[, 1] %in% commProSpe[, 6], ]
  KOList[[i]] <- KOMat
}

# whole species vector
wholeSpe = vector()
for (i in 1:8) {
  wholeSpe = union(wholeSpe, KOList[[i]][, 1])
}

# merge
ATPKO <- vector('list', length(wholeSpe))
names(ATPKO) <- wholeSpe
for (i in 1:length(wholeSpe)) {
  eachSpe <- lapply(KOList, function(x) {
    eachSpeKO <- x[x[, 1] %in% wholeSpe[i], ,drop = FALSE]
    # some species may lack certain subuints
    speKONum <- nrow(eachSpeKO)
    if (speKONum == 0) {
      uniKO <- NA
    } else {
      uniKO <- eachSpeKO[, 2]
    }

    return(unname(uniKO))
  })

  ATPKO[[i]] <- unlist(eachSpe)

}

# transfer KEGG ID to BioCyc ID
# BioCyc speID
uniSpe <- commProSpe[commProSpe[, 6] %in% wholeSpe, ]
uniSpe <- uniSpe[order(as.character(uniSpe[, 6])), ]
uniSpe <- uniSpe[rank(wholeSpe), ]
ATPKOKEGGSpe <- names(ATPKO)
names(ATPKO) <- as.character(uniSpe[, 2])

tmp1 <- ATPKO

# some may contain ''', like 'atpF'', and could not be identified by R ramote server
for (i in 1:length(ATPKO)) {
  print(paste('It is running ', i, '.', sep = ''))
  x <- ATPKO[[i]]
  KEGGIDVec <- paste(ATPKOKEGGSpe[i], x, sep = ':')
  KEGGIDVec <- paste(KEGGIDVec, collapse='+')
  KEGGsymTable <- webTable(paste('http://rest.kegg.jp/list/', KEGGIDVec, sep = ''), n = 2)
  KEGGsym <- KEGGsymTable[, 2]
  KEGGsym <- sapply(strsplit(KEGGsym, split = ';', fixed = TRUE), '[', 1)
  y <- character(length = length(x))
  y[which(is.na(x))] <- NA
  y[which(!is.na(x))] <- KEGGsym
  hasSym <- which(!grepl(' ', y))
  x[hasSym] <- y[hasSym]
  ATPKO[[i]] <- x
}

identical(sapply(tmp1, length), sapply(ATPKO, length))
identical(sapply(tmp1, is.na), sapply(tmp1, is.na))

# sperate ATPKO
hasDot <- sapply(ATPKO, function(x) {
  hasDotEach <- grepl('\'', x)
  if (sum(hasDotEach) > 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
})

save(ATPKO, ATPKOKEGGSpe, file = 'ATPKO.RData')

# ========================== script =======================
ATPKODot <- ATPKO[hasDot]
ATPKOKEGGSpeDOT <- ATPKOKEGGSpe[hasDot]
save(ATPKODot, ATPKOKEGGSpeDOT, file = 'ATPKODOT.RData')

ATPKO <- ATPKO[!hasDot]
ATPKOKEGGSpe <-ATPKOKEGGSpe[!hasDot]
save(ATPKO, ATPKOKEGGSpe, file = 'ATPKONODOT.RData')
# =========================================================


# BioCyc geneID. Cut the whole length with internal 4
cutMat <- CutSeqEqu(length(ATPKODot), 4)
for (j in 1:ncol(cutMat)) {
  ATPKOCycPart <- foreach (i = cutMat[1, j]:cutMat[2, j]) %dopar% {
    # may have NA
    print(paste('It is running ', i, ' with the name of ', names(ATPKODot)[i], '.', sep = ''))
    iniVal <- ATPKODot[[i]]
    iniVal[!is.na(iniVal)] <- sapply(iniVal[!is.na(iniVal)], KEGGID2CycID, speKEGGID = ATPKOKEGGSpeDOT[i], speCycID = names(ATPKODot)[i])
    return(iniVal)
  }
  names(ATPKOCycPart) <- names(ATPKODot)[cutMat[1, j]:cutMat[2, j]]
  save(ATPKOCycPart, file = paste('ATPKOCycPartDOT', cutMat[1, j], '_', cutMat[2, j], '.RData', sep = ''))
  Sys.sleep(30)
}

# BioCyc geneID. Cut the whole length with internal 4
cutMat <- CutSeqEqu(length(ATPKO), 4)
for (j in 1:ncol(cutMat)) {
  ATPKOCycPart <- foreach (i = cutMat[1, j]:cutMat[2, j]) %dopar% {
    # may have NA
    print(paste('It is running ', i, ' with the name of ', names(ATPKO)[i], '.', sep = ''))
    iniVal <- ATPKO[[i]]
    iniVal[!is.na(iniVal)] <- sapply(iniVal[!is.na(iniVal)], KEGGID2CycID, speKEGGID = ATPKOKEGGSpe[i], speCycID = names(ATPKO)[i])
    return(iniVal)
  }
  names(ATPKOCycPart) <- names(ATPKO)[cutMat[1, j]:cutMat[2, j]]
  save(ATPKOCycPart, file = paste('ATPKOCycNODOTPart', cutMat[1, j], '_', cutMat[2, j], '.RData', sep = ''))
  Sys.sleep(60)
}

list2list <- function(lsInput){
  # not support NA
  ## $alphaKO
  ## [1] "TU0-6636" "TU0-6635" "TU00243"

  ## $betaKO
  ## [1] "TU0-6636" "TU0-6635" "TU00243"

  ## $gammaKO
  ## [1] "TU0-6636" "TU0-6635" "TU00243"

  ## $detaKO
  ## [1] "TU0-6636" "TU0-6635" "TU00243"

  ## $epsilonKO
  ## [1] "TU0-42328" "TU0-6636"  "TU0-6635"  "TU00243"

  ## $cKO
  ## [1] "TU0-6636" "TU0-6635" "TU00243"

  ## $aKO
  ## [1] "TU0-6636" "TU0-6635" "TU00243"

  ## $bKO
  ## [1] "TU0-6636" "TU0-6635" "TU00243"


  uniEle <- unique(unlist(lsInput))
  uniList <- vector('list', length(uniEle))
  names(uniList) <- uniEle
  for (i in 1:length(uniEle)) {
    uniList[[i]] <- names(lsInput)[sapply(lsInput, function(x){uniEle[i] %in% x})]
  }
  return(uniList)
}

# transcription unit
ATPTU <- vector('list', length(ATPKO))
names(ATPTU) <- names(ATPKO)
for (i in 1:length(ATPKO)) {
  # get TU name
  eachTU <- lapply(ATPKO[[i]], function(x) {
    if (is.na(x)) {
      TU <- NA
    } else {
      eachGeneInfo <- getCycTUfGene(x, speID = names(ATPKO)[i])
      TU <- eachGeneInfo
      if (is.null(TU)) {
        TU <- NA
      } else {}
    }
    return(TU)
  })

  # range the list
  # select non NA
  hasNA <- sapply(eachTU, function(x) {
    if (sum(is.na(x)) > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  TUnona <- eachTU[which(!hasNA)]
  TUList <- list2list(TUnona)

  # select NA
  TUna <- eachTU[which(hasNA)]
  if (length(TUna) != 0) {
    TUList$noTU <- names(TUna)
  } else {}

  ATPTU[[i]] <- TUList
}

save(ATPTU, file = 'ATPCycPhyloTU.RData')

#################### process the repeat data ################
load('wholeTUList.RData')
delNum <- 284
nonDelNum <- which(!((1:length(wholeTUList)) %in% delNum))
wholeTUList <- wholeTUList[nonDelNum]

duNames <- names(wholeTUList)[which(duplicated(names(wholeTUList)))]
which(names(wholeTUList) %in% duNames)
wholeTUList[which(names(wholeTUList) %in% duNames)]


load('ATPKOCyc.RData')
delNum <- c(1240, 1293, 1428, 1857)
nonDelNum <- which(!((1:length(ATPKOCyc)) %in% delNum))
ATPKOCyc <- ATPKOCyc[nonDelNum]

duNames <- names(ATPKOCyc)[which(duplicated(names(ATPKOCyc)))]
which(names(ATPKOCyc) %in% duNames)
ATPKOCyc[which(names(ATPKOCyc) %in% duNames)]




############## calculate the loss/repeat ATP ############
ATPSubMat <- matrix(ncol = 8, nrow = length(ATPKO))
colnames(ATPSubMat) <- names(KOVec)
rownames(ATPSubMat) <- names(ATPKO)

for (i in 1:length(ATPKO)) {
  subGene <- ATPKO[[i]]
  subName <- names(subGene)
  ATPsubName <- names(KOVec)
  subNamePat <- paste('^', ATPsubName, 'KO', sep = '')
  subNum <- integer(8)
  for (j in 1:8) {
    sub <- subGene[grep(subNamePat[j], subName)]
    if (sum(is.na(sub)) > 0) {
      num <- 0
    } else {
      num <- length(sub)
    }
    subNum[j] <- num
  }
  ATPSubMat[i, ] <- subNum
}

write.csv(ATPSubMat, 'ATPSubMat.csv')


ATPSubMat2 <- apply(ATPSubMat, 1:2, function(x) {
  if (x > 0) {
    x <- 1
  } else {
    x <- 0
  }
  return(x)
})

require(pheatmap)
subComplexCol <- data.frame(Subunit = factor(c(rep(1, 5), rep(2, 3)), labels = c('F1', 'Fo')))
rownames(subComplexCol) <- names(KOVec)
pheatmap(as.matrix(dist(t(ATPSubMat), diag = TRUE, upper = TRUE)), annotation = subComplexCol, clustering_method = "average", cellwidth = 15, cellheight = 12)

########################## test code ########################

## test2 <- ATPKO[1:2]
## for (i in 1:length(test2)){
##   may have NA
##   iniVal <- test2[[i]]
##   iniVal[!is.na(iniVal)] <- names(KEGGID2CycID(iniVal[!is.na(iniVal)], names(test2)[i], n = 4))
##   test2[[i]] <- iniVal
## }


## test1 <- vector('list', 3)
## names(test1) <- c('ECOLI', 'RRUB269796', 'MTUB1304279-WGS')
## ecoKO <- c('EG10098', 'EG10101', 'EG10104', 'EG10105', 'EG10100', 'EG10102', 'EG10099', 'EG10103')
## names(ecoKO) <- names(ATPKO$ECOLI)
## rruKO <- c('GCN1-1247', 'GCN1-1249', 'GCN1-1248', 'GCN1-1246', 'GCN1-1250', 'GCN1-3300', 'GCN1-3301', 'GCN1-3298', 'GCN1-3299')
## names(rruKO) <- names(ATPKO$RRUB269796)
## mtuhKO <- c(NA, NA, 'GSRF-1250', NA, 'GSRF-1251', 'GSRF-1247', 'GSRF-1246', 'GSRF-1248')
## names(mtuhKO) <- names(ATPKO$`MTUB1304279-WGS`)
## test1[[1]] <- ecoKO
## test1[[2]] <- rruKO
## test1[[3]] <- mtuhKO



## eachTU <- lapply(test1[[2]], function(x) {
##     if (is.na(x)) {
##       TU <- NA
##     } else {
##       eachGeneInfo <- getCycTUfGene(x, speID = 'RRUB269796')
##       TU <- eachGeneInfo
##       if (is.null(TU)) {
##         TU <- NA
##       } else {}
##     }
##     return(TU)
##   })


## ATPTU <- vector('list', length(test1))
## names(ATPTU) <- names(test1)
## for (i in 1:length(test1)) {
##   get TU name
##   eachTU <- lapply(test1[[i]], function(x) {
##     if (is.na(x)) {
##       TU <- NA
##     } else {
##       eachGeneInfo <- getCycTUfGene(x, speID = names(test1)[i])
##       TU <- eachGeneInfo
##       if (is.null(TU)) {
##         TU <- NA
##       } else {}
##     }
##     return(TU)
##   })

##   range the list
##   select non NA
##   hasNA <- sapply(eachTU, function(x) {
##     if (sum(is.na(x)) > 0) {
##       return(TRUE)
##     } else {
##       return(FALSE)
##     }
##   })

##   TUnona <- eachTU[which(!hasNA)]
##   TUList <- list2list(TUnona)

##   select NA
##   TUna <- eachTU[which(hasNA)]
##   if (length(TUna) != 0) {
##     TUList$noTU <- names(TUna)
##   } else {}

##   ATPTU[[i]] <- TUList
## }
