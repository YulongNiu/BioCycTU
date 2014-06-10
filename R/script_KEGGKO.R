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
ATPCycList$MBOU1201294['A'] <- 'GLGD-172'
ATPCycList$MBOU1201294['F'] <- 'GLGD-171'
ATPCycList$MBOU1201294['K'] <- 'GLGD-168'
ATPCycList$`CTRA1071771-WGS`['I'] <- 'GSK2-309'
# delete TOSH751945
ATPCycList <- ATPCycList[!(names(ATPCycList) %in% 'TOSH751945')]
# duplicated names
ATPCycList <- ATPCycList[!duplicated(names(ATPCycList))]
save(ATPCycList, file = 'ATPCycList.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############### TU ###################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('ATPCycList.RData')
ATPTU <- vector('list', length(ATPCycList))
names(ATPTU) <- names(ATPCycList)
for (i in 66:length(ATPCycList)) {
  print(paste('It is running ', i, '.', sep = ''))
  # get TU name
  eachTU <- lapply(ATPCycList[[i]], function(x) {
    if (is.na(x)) {
      TU <- NA
    } else {
      eachGeneInfo <- getCycTUfGene(x, speID = names(ATPCycList)[i])
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
  Sys.sleep(30)
}
ATPTU <- ATPTU[!duplicated(names(ATPTU))]
save(ATPTU, file = 'ATPTU.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('ATPTU.RData')
load('ATPCycList.RData')
for (i in 508:length(ATPTU)) {
  print(paste('It is running ', i, ' with the name ', names(ATPTU[i]), sep = ''))
  TUListOrder <- ATPTU[[i]]
  speID <- names(ATPTU[i])
  speGenes <- ATPCycList[[which(names(ATPCycList) %in% speID)]]

  for (j in 1:length(ATPTU[[i]])) {
    if (names(ATPTU[[i]][j]) != 'noTU') {
      if (length(ATPTU[[i]][[j]]) > 1) {
        cycGenesOrder <- getCycTUInfo(names(ATPTU[[i]])[j], speID)
        cycGenesOrder <- cycGenesOrder$geneIDs
        TUKO <- ATPTU[[i]][[j]]
        TUgenes <- speGenes[names(speGenes) %in% TUKO]
        TUgenesOrder <- cycGenesOrder[cycGenesOrder %in% TUgenes]
        TUgenesNames <- sapply(TUgenesOrder, function(x) {
          y <- names(speGenes[speGenes %in% x])
          return(y)
        })
        names(TUgenesOrder) <- TUgenesNames
      }
      else if (length(ATPTU[[i]][[j]]) == 1){
        # only have one element
        TUgenesOrder <- speGenes[names(speGenes) %in% ATPTU[[i]][[j]]]
      }
    } else {
      TUgenesOrder = speGenes[names(speGenes) %in% ATPTU[[i]][[j]]]
    }

    TUListOrder[[j]] <- TUgenesOrder
  }
  TUListOrder <- list(TUListOrder)
  names(TUListOrder) <- names(ATPTU[i])

  save(TUListOrder, file = paste('TUListOrder_', i, '.RData', sep = ''))

}

TUOrder <- dir(pattern="^TUListOrder_.*")
TUOrderList <- foreach(i = 1:length(TUOrder), .combine = append) %dopar% {
  load(TUOrder[i])
  return(TUListOrder)
}

save(TUOrderList, file = 'TUOrderList.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################## continuous TU #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('TUOrderList.RData')
TUOrderListCon <- vector('list', length = length(TUOrderList))
names(TUOrderListCon) <- names(TUOrderList)

contiMat <- function(vec) {
  # example contiMat(c(1:3, 5:7, 9))
  diffs <- c(1, diff(vec))
  start_indexes <- c(1, which(abs(diffs) > 1))
  end_indexes <- c(start_indexes - 1, length(vec))
  contiMat <- cbind(vec[start_indexes], vec[end_indexes])
  return(contiMat)
}

TUListEachCon <- function(x) {
  # x is an list of TU
  TUNames <- names(x)
  if ('noTU' %in% TUNames) {
    noTUList <- list(noTU = x$noTU)
    TUNames <- TUNames[!(TUNames %in% 'noTU')]
  } else {
    noTUList <- list(noTU = NULL)
  }

  # split names by '-'
  TUNamesNum <- sapply(strsplit(TUNames, split = '-', fixed = TRUE), '[', 2)
  names(TUNames) <- TUNamesNum
  TUNamesNum <- as.numeric(TUNamesNum)

  # order TU numbers
  TUNamesMat <- contiMat(sort(TUNamesNum))

  TUListCon <- vector('list', length = nrow(TUNamesMat))

  for (i in 1:nrow(TUNamesMat)) {
    if (diff(c(TUNamesMat[i, 1], TUNamesMat[i, 2])) == 0) {
      selectTUNames <- names(TUNames) %in% as.character(TUNamesMat[i, 1])
      TUListCon[i] <- x[which(selectTUNames)]
      names(TUListCon)[i] <- TUNames[which(selectTUNames)]
    } else {
      selectTUNames <- names(TUNames) %in% as.character(TUNamesMat[i, 1] : TUNamesMat[i, 2])
      TUListConSel <- unlist(x[which(selectTUNames)])
      namesTUListConSel <- names(TUListConSel)
      namesTUListConSel <- sapply(strsplit(namesTUListConSel, split = '.', fixed = TRUE), '[', 2)
      names(TUListConSel) <- namesTUListConSel
      TUListCon[[i]] <- sort(TUListConSel)
      names(TUListCon)[i] <- paste(TUNames[which(selectTUNames)], collapse = '|')
    }
  }

  if (!is.null(noTUList$noTU)) {
    TUListCon <- append(TUListCon, noTUList)
  } else {}

  return(TUListCon)
}

TUOrderListCon <- sapply(TUOrderList, TUListEachCon)
save(TUOrderListCon, file = 'TUOrderListCon.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

########### precess the TUOrder data, add "KO" in each subunits' name ##########
load('TUOrderList.RData')
load('TUOrderListCon.RData')

getcontent <- function(s,g) {
  substring(s,g,g+attr(g,'match.length')-1)
}

TUOrderList <- lapply(TUOrderList, function(x) {
  lapply(x, function(y) {
    # select the character ones
    chaName <- regexpr('[a-zA-Z]+', names(y))
    chaName <- getcontent(names(y), chaName)
    numName <- regexpr('\\d+', names(y))
    numName <- getcontent(names(y), numName)
    names(y) <- paste(chaName, 'KO', numName, sep = '')
    return(y)
  })
})

TUOrderListCon <- lapply(TUOrderListCon, function(x) {
  lapply(x, function(y) {
    # select the character ones
    chaName <- regexpr('[a-zA-Z]+', names(y))
    chaName <- getcontent(names(y), chaName)
    numName <- regexpr('\\d+', names(y))
    numName <- getcontent(names(y), numName)
    names(y) <- paste(chaName, 'KO', numName, sep = '')
    return(y)
  })
})

save(TUOrderList, file = 'TUOrderList.RData')
save(TUOrderListCon, file = 'TUOrderListCon.RData')

################## process TU order ################
load('TUOrderList.RData')
load('TUOrderListCon.RData')

stanCut <- c('D', 'B', 'A', 'F', 'C', 'E', 'K', 'I', 'H')

preTUCut <- function(preTU, stanSeq, noTU = FALSE){
  # pre-process the input TU sequence
  # if noTU is TRUE, we prefer to order the noTU elements according to the stanSeq.
  greStanSeq <- paste('^', stanSeq, 'KO', sep = '')
  TU <- vector(length = length(preTU))

  getcontent <- function(s,g) {
    substring(s,g,g+attr(g,'match.length')-1)
  }

  for (i in 1:length(preTU)){
    for (j in 1:length(greStanSeq)) {
      grepNum <- gregexpr(greStanSeq[j], preTU[i])[[1]]
      if (grepNum > 0) {
        grepSeq <- getcontent(preTU[i], grepNum)
        grepSeq <- substr(grepSeq, 1, nchar(grepSeq) - 2)
        TU[i] <- grepSeq
      } else {}
    }
  }

  #========= remove repeat, beacuse I suppose the repeat doesn't has matter to cut point ===========
  TU <- TU[!duplicated(TU)]
  if (noTU) {
    stanCutTU <- stanCut[stanCut %in% TU]
    # order noTU elements
    TU <- TU[order(TU)]
    TU <- TU[rank(stanCutTU)]
  } else {}

  return(TU)

}

##' Get the cut point by TU.
##'
##' Get the cut point by a given transcription unit
##' @title Get the cut point.
##' @param TU The input TU vector.
##' @param stanSeq The standard TU sequence
##' @return A vector of the cut point
##' @examples
##' @author Yulong Niu \email{niuylscu@@gmail.com}
cutTU <- function(TU, stanSeq) {

  TUNum <- sapply(TU, function(x){
    num <- which(stanSeq %in% x)
    return(num)
  })

  contiMat <- function(vec) {
    # example contiMat(c(1:3, 5:7, 9))
    diffs <- c(1, diff(vec))
    start_indexes <- c(1, which(abs(diffs) > 1))
    end_indexes <- c(start_indexes - 1, length(vec))
    contiMat <- cbind(vec[start_indexes], vec[end_indexes])
    return(contiMat)
  }

  cutTUmat <- contiMat(TUNum)

  # +1 to the bigger cut point in each row
  cutTUmat <- t(apply(cutTUmat, 1, function(x){
    x[which.max(x)] <- x[which.max(x)] + 1
    return(x)
  }))

  return(cutTUmat)

}



load('TUOrderList.RData')
library(foreach)
library(doMC)
registerDoMC(4)

CutListVec <- function(orderList, stanSeq, cutIn = FALSE) {
  TUCutList <- lapply(orderList, function(x) {
    eachTU <- vector('list', length = length(x))
    for (i in 1:length(x)) {
      names(eachTU) <- names(x)
      if (names(x[i]) == 'noTU') {
        eachTUPreCut <- preTUCut(names(x[[i]]), stanSeq, noTU = TRUE)
      } else {
        eachTUPreCut <- preTUCut(names(x[[i]]), stanSeq)
      }
      eachTU[[i]] <- cutTU(eachTUPreCut, stanSeq)
    }
    return(eachTU)
  })

  if (!cutIn) {
    TUCutMat <- foreach(i = 1:length(TUCutList), .combine = rbind) %dopar% {
      foreach(j = 1:length(TUCutList[[i]]), .combine = rbind) %dopar% {
        t(apply(TUCutList[[i]][[j]], 1, sort))
      }
    }
  } else {
    # summary the cut in number in one TU
    TUCutMat <- foreach(i = 1:length(TUCutList), .combine = rbind) %dopar% {
      foreach(j = 1:length(TUCutList[[i]]), .combine = rbind) %dopar% {
        if (nrow(TUCutList[[i]][[j]]) > 1) {
          t(apply(TUCutList[[i]][[j]], 1, sort))
        } else {}
      }
    }
  }

  TUCutVec <- apply(TUCutMat, 1, paste, collapse = '-')

  TUCutTable <- table(TUCutVec)
  fullName <- apply(combn(length(stanCut)+1, 2), 2, paste, collapse = '-')
  fullIni <- numeric(length = length(fullName))
  names(fullIni) <- fullName
  fullIni[which(names(fullIni) %in% sort(names(TUCutTable)))] <- TUCutTable

  ## cutNum <- sapply(TUCutVec, function(x) {
  ##   eachTUCut <- sum(sapply(x, nrow))
  ## })
  ## TUNum <- sapply(orderList, length)

  ## return(fullIni)
  return(fullIni)
}

TUcutNum <- CutListVec(TUOrderListCon, stanCut)
