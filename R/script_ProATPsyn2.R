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
load('ATPKONODOT.RData')
library(foreach)
library(doMC)
library(XML)
registerDoMC(1)

# some species name changed
names(ATPKO)[names(ATPKO) %in% 'BANT198094-WGS'] <- 'ANTHRA'
names(ATPKO)[names(ATPKO) %in% 'CCRE190650-WGS'] <- 'CAULO'
names(ATPKO)[names(ATPKO) %in% 'MTUB83332-WGS'] <- 'MTBRV'


# BioCyc geneID. Cut the whole length with internal 4
cutMat <- CutSeqEqu(length(ATPKO), 4)
for (j in 452:ncol(cutMat)) {
  ATPKOCycPart <- foreach (i = cutMat[1, j]:cutMat[2, j]) %dopar% {
    # may have NA
    print(paste('It is running ', i, ' with the name of ', names(ATPKO)[i], '.', sep = ''))
    iniVal <- ATPKO[[i]]
    iniVal[!is.na(iniVal)] <- sapply(iniVal[!is.na(iniVal)], KEGGID2CycID, speKEGGID = ATPKOKEGGSpe[i], speCycID = names(ATPKO)[i])
    return(iniVal)
    Sys.sleep(30)
  }
  names(ATPKOCycPart) <- names(ATPKO)[cutMat[1, j]:cutMat[2, j]]
  save(ATPKOCycPart, file = paste('ATPKONODOTCycPart', cutMat[1, j], '_', cutMat[2, j], '.RData', sep = ''))
  if (j %% 20 == 0) {
    Sys.sleep(300)
  } else {}
}
