##' Cut the index according to the input pattern
##'
##' Cut the index of a vector based on the input pattern. There is not need to input the index. The sum of 'cutSeq' should be equal to the length of vector, which is used to cut.
##' @title Cut index
##' @param cutSeq The input cut pattern. 0 is automatically removed.
##' @return A matrix of index with two rows.
##' @examples CutSeq(c(1, 4, 8))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
CutSeq <- function(cutSeq){

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
##' Cut the index with equal internal.
##'
##' The index is cut with equal internals.
##' @title Cut index with the equal internal.
##' @param vecLen The length of vector used to cut.
##' @param equNum Internal number
##' @return A index matrix.
##' @examples
##' CutSeqEqu(10, 2)
##' CutSeqEqu(10, 3)
##' # the internal is bigger the length of input index
##' CutSeqEqu(5, 6)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
CutSeqEqu <- function(vecLen, equNum){

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



##' Rearrange KO by species
##'
##' It rearrange the given KO vectors according to the species. If the KO gene is missing in one species, it is marked as NA.
##' @title Rearrange KO by species ID
##' @param KOvec Vector of KO
##' @param KOnames names of KO
##' @param n The number of CPUs or processors, and the default value is 4.
##' @return
##' @examples
##' # Bacterial F1Fo ATP synthase KO
##' FATPKO <- c('K02111', 'K02112', 'K02115', 'K02113', 'K02114', 'K02110', 'K02108', 'K02109')
##' FATPKONames <- c('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'c', 'a', 'b')
##' FATPKObySpe <- KEGGSpeKO(FATPKO, FATPKONames)
##' # Bacterial A1Ao ATP synthase KO
##' AATPKO <- c('K02117', 'K02118', 'K02119', 'K02120', 'K02121', 'K02122', 'K02107', 'K02123', 'K02124')
##' AATPKONames <- c('A', 'B', 'C', 'D', 'E', 'F', 'H', 'I', 'K')
##' AATPKObySpe <- KEGGSpeKO(AATPKO, AATPKONames)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom KEGGBioCycAPI getKEGGKO
##' @importFrom doMC registerDoMC
##' @importFrom foreach foreach
##'
KEGGSpeKO <- function(KOvec, KOnames = NULL, n = 4) {

  require(KEGGBioCycAPI)
  require(foreach)
  require(doMC)
  registerDoMC(n)

  if (is.null(KOvec)) {
    # without names, so use the KOvec as the names
    KOnames <- KOvec
  } else {
    if (!identical(length(KOvec), length(KOnames))) {
      stop('The input "KOvec" and "KOnames" should be in the same length.')
    } else {}
  }

  # parallel getting KO list
  KOListPara <- foreach(i = 1:length(KOvec), .inorder = FALSE) %dopar% {
    KOgenes <- getKEGGKO(KOvec[i])
    KOMat <- unlist(strsplit(KOgenes, split = ':', fixed = TRUE))
    KOMat <- matrix(KOMat, ncol = 2, byrow = TRUE)
    KOnm <- KOnames[i]
    eachKO <- list(KOMat = KOMat, KOnm = KOnm)
    return(eachKO)
  }

  # KOList with KEGG species ID
  KOList <- sapply(KOListPara, '[[', 1)
  names(KOList) <- sapply(KOListPara, '[[', 2)

  # whole species List
  wholeSpe = vector()
  for (i in 1:8) {
    wholeSpe = union(wholeSpe, KOList[[i]][, 1])
  }

  speKOList <- vector('list', length(wholeSpe))
  names(speKOList) <- wholeSpe
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

    speKOList[[i]] <- unlist(eachSpe)
  }

  return(speKOList)
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
