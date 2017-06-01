# devtools::install_github('yufree/sva-devel')
# devtools::install_github('yufree/enviGCMS')
library(enviGCMS)
library(BioMark)

data("SpikePos")
data("SpikeNeg")
# function re-write for BioMark package
svacorBM <-
        function(SpikePos,
                 lv = NULL) {
                data <- t(SpikePos$data)
                
                mz <- SpikePos$annotation$mz
                rt <- SpikePos$annotation$rt
                mod <- stats::model.matrix(~ lv)
                mod0 <- as.matrix(c(rep(1, ncol(data))))
                svafit <- sva::sva(data, mod)
                if (svafit$n.sv == 0) {
                        svaX <- stats::model.matrix(~ lv)
                        lmfit <- limma::lmFit(data, svaX)
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal
                        rownames(signal) <-
                                rownames(error) <- rownames(data)
                        colnames(signal) <-
                                colnames(error) <- colnames(data)
                        pValues = sva::f.pvalue(data, mod, mod0)
                        qValues = qvalue::qvalue(pValues)
                        qValues = qValues$qvalues
                        
                        li <- list(data,
                                   signal,
                                   error,
                                   pValues,
                                   qValues,
                                   mz,
                                   rt)
                        names(li) <-
                                c(
                                        'data',
                                        'signal',
                                        'error',
                                        'p-values',
                                        'q-values',
                                        'mz',
                                        'rt'
                                )
                }
                else{
                        message('Data is correcting ...')
                        svaX <- stats::model.matrix(~ lv + svafit$sv)
                        lmfit <- limma::lmFit(data, svaX)
                        batch <-
                                lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + svafit$n.sv)] %*% t(svaX[, (nlevels(lv) +
                                                                                                                  1):(nlevels(lv) + svafit$n.sv)])
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal - batch
                        datacor <- signal + error
                        svaX2 <- stats::model.matrix(~ lv)
                        lmfit2 <- limma::lmFit(data, svaX2)
                        signal2 <-
                                lmfit2$coef[, 1:nlevels(lv)] %*% t(svaX2[, 1:nlevels(lv)])
                        error2 <- data - signal2
                        rownames(signal2) <-
                                rownames(error2) <-
                                rownames(datacor) <-
                                rownames(signal) <-
                                rownames(batch) <-
                                rownames(error) <- rownames(data)
                        colnames(signal2) <-
                                colnames(error2) <-
                                colnames(datacor) <-
                                colnames(signal) <-
                                colnames(batch) <-
                                colnames(error) <- colnames(data)
                        
                        modSv = cbind(mod, svafit$sv)
                        mod0Sv = cbind(mod0, svafit$sv)
                        pValuesSv = sva::f.pvalue(data, modSv, mod0Sv)
                        qValuesSv = qvalue::qvalue(pValuesSv)
                        qValuesSv = qValuesSv$qvalues
                        
                        pValues = sva::f.pvalue(data, mod, mod0)
                        qValues = qvalue::qvalue(pValues)
                        qValues = qValues$qvalues
                        
                        li <-
                                list(
                                        data,
                                        datacor,
                                        signal,
                                        batch,
                                        error,
                                        signal2,
                                        error2,
                                        pValues,
                                        qValues,
                                        pValuesSv,
                                        qValuesSv,
                                        svafit$pprob.gam,
                                        svafit$pprob.b,
                                        mz,
                                        rt
                                )
                        names(li) <-
                                c(
                                        'data',
                                        'dataCorrected',
                                        'signal',
                                        'batch',
                                        'error',
                                        'signal2',
                                        'error2',
                                        'p-values',
                                        'q-values',
                                        'p-valuesCorrected',
                                        'q-valuesCorrected',
                                        'PosteriorProbabilitiesSurrogate',
                                        'PosteriorProbabilitiesMod',
                                        'mz',
                                        'rt'
                                )
                        message('Done!')
                }
                return(li)
        }

# set the group
lv <- as.factor(c(rep('C',10),rep('1',10),rep('2',10),rep('3',10)))
lv <- relevel(lv, "C")

# SVA correction on positive and negative mode

## positive mode pca(figure 1)

li <- svacorBM(SpikePos = SpikePos, lv = lv)

pdf('pcap.pdf',width = 10.5, height = 7)
enviGCMS::svapca(li,lv = as.character(lv))
dev.off()

## negative mode pca(figure S1)

lin <- svacorBM(SpikePos = SpikeNeg, lv = lv)

pdf('pcan.pdf',width = 10.5, height = 7)
enviGCMS::svapca(lin,lv = as.character(lv))
dev.off()

# heatmap linear decomposition

## positive mode raw data(figure S2)

pdf('posheat.pdf', width = 8, height = 6)
dfanova <- enviGCMS::svaplot(li,pqvalues = 'anova', lv = lv)
dev.off()

## positive mode corrected data(figure 2)

pdf('posheatsv.pdf', width = 12, height = 6)
dfsv <- enviGCMS::svaplot(li,pqvalues = 'sv', lv = lv)
dev.off()

## negative mode raw data(figure S3)

pdf('posheatn.pdf', width = 8, height = 6)
dfanovan <- enviGCMS::svaplot(lin,pqvalues = 'anova', lv = lv)
dev.off()

## negative mode raw data(figure S4)

pdf('posheatsvn.pdf', width = 12, height = 6)
dfsvn <- enviGCMS::svaplot(lin,pqvalues = 'sv', lv = lv)
dev.off()

## positive mode corrected data for pre-defined biomarkers

# get the index

ind <- SpikePos$annotation$found.in.standards != 0
indn <- SpikeNeg$annotation$found.in.standards != 0

## postive mode for raw data

pdf('posheatsub.pdf', width = 8, height = 6)
dfanovasub <- enviGCMS::svaplot(li,pqvalues = 'anova', lv = lv, index = ind)
dev.off()

## postive mode for corrected data

pdf('posheatsvsub.pdf', width = 12, height = 6)
dfsvsub <- enviGCMS::svaplot(li,pqvalues = 'sv', lv = lv,index = ind)
dev.off()

## negative mode for raw data

pdf('posheatsubn.pdf', width = 8, height = 6)
dfanovansub <- enviGCMS::svaplot(lin,pqvalues = 'anova', lv = lv, index = indn)
dev.off()

## negative mode for corrected data

pdf('posheatsvsubn.pdf', width = 12, height = 6)
dfsvnsub <- enviGCMS::svaplot(lin,pqvalues = 'sv', lv = lv,index = indn)
dev.off()

# plot the influnces on DoE and batch effects

## positive mode (figure 3)

pdf('relative.pdf', height = 5, width = 8)
enviGCMS::svabatch(li,dfsv,dfanova)
dev.off()

## negative mode (figure S5)

pdf('relativen.pdf', height = 5, width = 8)
enviGCMS::svabatch(lin,dfsvn,dfanovan)
dev.off()

# quantitative analysis in positive mode

## get the pre-defined biomarks peaks found statistically significant in raw and corrected data from raw data

dataraw <- t(SpikePos$data)[SpikePos$annotation$found.in.standards != 0&dfsv$pqvalues&dfanova$pqvalues,]

## get the pre-defined biomarks peaks found statistically significant in raw and corrected data from corrected data

datacor <- li$dataCorrected[SpikePos$annotation$found.in.standards != 0&dfsv$pqvalues&dfanova$pqvalues,] 

## Combine them
datax <- cbind(dataraw,datacor)

## Get their m/z and RT

mzs <- SpikePos$annotation$mz[SpikePos$annotation$found.in.standards != 0&dfsv$pqvalues&dfanova$pqvalues]

rts <- SpikePos$annotation$rt[SpikePos$annotation$found.in.standards != 0&dfsv$pqvalues&dfanova$pqvalues]

## Get the influnces from surrogate variable and primary variable

p0 <- li$PosteriorProbabilitiesSurrogate[SpikePos$annotation$found.in.standards != 0&dfsv$pqvalues&dfanova$pqvalues]
p1 <- li$PosteriorProbabilitiesMod[SpikePos$annotation$found.in.standards != 0&dfsv$pqvalues&dfanova$pqvalues]

## Get the distance ratio for those peaks
dsm <- matrix(nrow = nrow(datax), ncol = 6)
for(i in 1:nrow(datax)){
        datam <- matrix(as.numeric(datax[i,]),10,8)
        mean <- apply(datam,2,mean)
        mean1 <- mean[2:4]
        mean2 <- mean[6:8]
        ds <- as.numeric(dist(mean1))[order(as.numeric(dist(mean1)))]/min(as.numeric(dist(mean1))[order(as.numeric(dist(mean1)))])
        dsc <- as.numeric(dist(mean2))[order(as.numeric(dist(mean2)))]/min(as.numeric(dist(mean2))[order(as.numeric(dist(mean2)))])
        dsm[i,] <- c(ds,dsc)
}

## Combine them with annotation
dsr <- cbind.data.frame(dsm,mzs,rts,p0,p1)
dsr$mz <- round(dsr$mzs,1)
dsr$rt <- round(dsr$rts,-2)
posm <- pos.markers[complete.cases(pos.markers),]
posm$mz <- round(posm$mz,1)
posm$rt <- round(posm$rt,-2)
dsc <- merge(dsr,posm,by = c('mz','rt'))

## Get the mean distance ratio
### raw data
mean(dsc$`2`)
mean(dsc$`3`)
### corrected data
mean(dsc$`5`)
mean(dsc$`6`)

# quantitative analysis in negative mode

## get the pre-defined biomarks peaks found statistically significant in raw and corrected data from raw data

datarawn <- t(SpikeNeg$data)[SpikeNeg$annotation$found.in.standards != 0&dfsvn$pqvalues&dfanovan$pqvalues,]

## get the pre-defined biomarks peaks found statistically significant in raw and corrected data from corrected data

datacorn <- lin$dataCorrected[SpikeNeg$annotation$found.in.standards != 0&dfsvn$pqvalues&dfanovan$pqvalues,] 

## Combine them
dataxn <- cbind(datarawn,datacorn)

## Get their m/z and RT
mzsn <- SpikeNeg$annotation$mz[SpikeNeg$annotation$found.in.standards != 0&dfsvn$pqvalues&dfanovan$pqvalues]

rtsn <- SpikeNeg$annotation$rt[SpikeNeg$annotation$found.in.standards != 0&dfsvn$pqvalues&dfanovan$pqvalues]

## Get the influnces from surrogate variable and primary variable
p0 <- lin$PosteriorProbabilitiesSurrogate[SpikeNeg$annotation$found.in.standards != 0&dfsvn$pqvalues&dfanovan$pqvalues]
p1 <- lin$PosteriorProbabilitiesMod[SpikeNeg$annotation$found.in.standards != 0&dfsvn$pqvalues&dfanovan$pqvalues]

## Get the distance ratio for those peaks
dsmn <- matrix(nrow = nrow(dataxn), ncol = 6)
for(i in 1:nrow(dataxn)){
        datam <- matrix(as.numeric(dataxn[i,]),10,8)
        mean <- apply(datam,2,mean)
        mean1 <- mean[2:4]
        mean2 <- mean[6:8]
        ds <- as.numeric(dist(mean1))[order(as.numeric(dist(mean1)))]/min(as.numeric(dist(mean1))[order(as.numeric(dist(mean1)))])
        dsc <- as.numeric(dist(mean2))[order(as.numeric(dist(mean2)))]/min(as.numeric(dist(mean2))[order(as.numeric(dist(mean2)))])
        dsmn[i,] <- c(ds,dsc)
}

## Combine them with annotation
dsrn <- cbind.data.frame(dsmn,mzsn,rtsn,p0,p1)
dsrn$mz <- round(dsrn$mzsn,1)
dsrn$rt <- round(dsrn$rtsn,-2)
negm <- neg.markers[complete.cases(neg.markers),]
negm$mz <- round(negm$mz,1)
negm$rt <- round(negm$rt,-2)
dscn <- merge(dsrn,negm,by = c('mz','rt'))

## Get the mean distance ratio
### raw data
mean(dscn$`2`)
mean(dscn$`3`)
### corrected data
mean(dscn$`5`)
mean(dscn$`6`)

