devtools::install_github('yufree/sva-devel')
suppressWarnings(suppressPackageStartupMessages(library(xcms)))
library(RColorBrewer)
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
library(CAMERA)
library(qvalue)
#' Get xcmsset object in one step with optimized methods.
#' path the path to your data
#' index the index of the files
#' BPPARAM used for BiocParallel package
#' pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
#' ... arguments for xcmsSet function
#' the parameters are extracted from the papers. If you use name other than the name above, you will use the default setting of XCMS. Also I suggest IPO packages or apLCMS packages to get reasonable data for your own instrumental. If you want to summit the results to a paper, remember to include those parameters.
#' return a xcmsset object for that path or selected samples
#' references Patti, G. J.; Tautenhahn, R.; Siuzdak, G. Nat. Protocols 2012, 7 (3), 508–516.
getdata <-
        function(path,
                 index = F,
                 BPPARAM = BiocParallel::SnowParam(workers = 12),
                 pmethod = 'hplcorbitrap',
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
                if (index) {
                        cdffiles <- cdffiles[index]
                }
                if (pmethod == 'hplcorbitrap') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(10, 60),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcorbitrap') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(5, 20),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        xset <-
                                xcms::group(xset, bw = 2, mzwid = 0.015)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected data
                        xset2 <-
                                xcms::group(xset2, bw = 2, mzwid = 0.015)
                        xset3 <-
                                xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                } else if (pmethod == 'hplcqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.025)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.025)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'hplchqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 2,
                                                    mzwid = 0.025)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 2,
                                                    mzwid = 0.025)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplchqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 2,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 2,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else{
                        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM, ...)
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <- xcms::group(xset2)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                }
                return(xset3)
        }
#' svacor is used to correct the data by surrogate variable analysis and the improved algorithm would be loaded from yufree/sva-devel repo
#' the impute data should be a xcmsSet object
#' lv means the group infomation, which could be extract from the xcmsSet
#' annotation means whether use CAMERA package to annotate the data
#' polarity and the nSlaves are the parameters for CAMERA package
#' the output of this function is a list depending on the correction process
#' When no surrogate variables found, list would be the linear decomposation of raw data.
#' When surrogate variables found, list would contain the linear decomposation of raw data and corrected data
#' Results also contained one way anova results(both p-values and q-values)according to the experiment design and local FDRs for each peaks

svacor <-
        function(xset,
                 lv = NULL,
                 annotation = F,
                 polarity = "positive",
                 nSlaves = 12) {
                data <- groupval(xset, "maxint", value = 'into')
                if (is.null(lv)) {
                        lv <- xset@phenoData[, 1]
                }
                mz <- xset@groups[, 1]
                rt <- xset@groups[, 4]
                mod <- model.matrix( ~ lv)
                mod0 <- as.matrix(c(rep(1, ncol(data))))
                svafit <- sva(data, mod)
                if (svafit$n.sv == 0) {
                        svaX <- model.matrix( ~ lv)
                        lmfit <- lmFit(data, svaX)
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal
                        rownames(signal) <-
                                rownames(error) <- rownames(data)
                        colnames(signal) <-
                                colnames(error) <- colnames(data)
                        pValues = f.pvalue(data, mod, mod0)
                        qValues = qvalue(pValues)
                        qValues = qValues$qvalues
                        if (annotation) {
                                dreport <-
                                        annotateDiffreport(
                                                xset,
                                                metlin = T,
                                                polarity = polarity,
                                                nSlaves = nSlaves
                                        )
                                dreport <-
                                        dreport[order(as.numeric(rownames(dreport))), ]
                                li <-
                                        list(data,
                                             signal,
                                             error,
                                             pValues,
                                             qValues,
                                             dreport,
                                             mz,
                                             rt)
                                names(li) <-
                                        c(
                                                'data',
                                                'signal',
                                                'error',
                                                'p-values',
                                                'q-values',
                                                'diffreport',
                                                'mz',
                                                'rt'
                                        )
                        } else{
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
                }
                else{
                        message('Data is correcting ...')
                        svaX <- model.matrix( ~ lv + svafit$sv)
                        lmfit <- lmFit(data, svaX)
                        batch <-
                                lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + svafit$n.sv)] %*% t(svaX[, (nlevels(lv) +
                                                                                                                  1):(nlevels(lv) + svafit$n.sv)])
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal - batch
                        datacor <- signal + error
                        svaX2 <- model.matrix( ~ lv)
                        lmfit2 <- lmFit(data, svaX2)
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
                        pValuesSv = f.pvalue(data, modSv, mod0Sv)
                        qValuesSv = qvalue(pValuesSv)
                        qValuesSv = qValuesSv$qvalues
                        
                        pValues = f.pvalue(data, mod, mod0)
                        qValues = qvalue(pValues)
                        qValues = qValues$qvalues
                        if (annotation) {
                                dreport <-
                                        annotateDiffreport(
                                                xset,
                                                metlin = T,
                                                polarity = polarity,
                                                nSlaves = nSlaves
                                        )
                                dreport <-
                                        dreport[order(as.numeric(rownames(dreport))), ]
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
                                                dreport,
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
                                                'diffreport',
                                                'PosteriorProbabilitiesSurrogate',
                                                'PosteriorProbabilitiesMod',
                                                'mz',
                                                'rt'
                                        )
                        }
                        else{
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
                        }
                        message('Done!')
                }
                return(li)
        }
#' svapca is used to show the decomposation of the raw data or correct data
#' the input data could be the list from svacor fuction
svapca <- function(list,
                   center = T,
                   scale = T,
                   lv = NULL) {
        data <- list$data
        Signal <- list$signal
        Batch <- list$batch
        error <- list$error
        datacor <- list$dataCorrected
        if(is.null(lv)){
                pch = colnames(data)
        }else{
                pch = lv
        }
        
        par(mfrow = c(2, 5), mar = c(4, 4, 2.6, 1))
        
        pcao <- prcomp(t(data), center = center, scale = scale)
        pcaoVars = signif(((pcao$sdev) ^ 2) / (sum((pcao$sdev) ^ 2)), 3) *
                100
        plot(pcao, type = "l", main = "PCA")
        
        pca <- prcomp(t(Signal), center = TRUE, scale = TRUE)
        pcaVars = signif(((pca$sdev) ^ 2) / (sum((pca$sdev) ^ 2)), 3) *
                100
        plot(pca, type = "l", main = "PCA-signal")
        
        pcab <- prcomp(t(Batch), center = center, scale = scale)
        pcabVars = signif(((pcab$sdev) ^ 2) / (sum((pcab$sdev) ^ 2)), 3) *
                100
        plot(pcab, type = "l", main = "PCA-batch")
        
        pcae <- prcomp(t(datacor), center = center, scale = scale)
        pcaeVars = signif(((pcae$sdev) ^ 2) / (sum((pcae$sdev) ^ 2)), 3) *
                100
        plot(pcae, type = "l", main = "PCA-error")
        
        pcac <- prcomp(t(datacor), center = center, scale = scale)
        pcacVars = signif(((pcac$sdev) ^ 2) / (sum((pcac$sdev) ^ 2)), 3) *
                100
        plot(pcac, type = "l", main = "PCA-corrected")
        
        plot(
                pcao$x[, 1],
                pcao$x[, 2],
                xlab = paste("PC1:", pcaoVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaoVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA"
        )
        
        plot(
                pca$x[, 1],
                pca$x[, 2],
                xlab = paste("PC1:", pcaVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-signal"
        )
        
        plot(
                pcab$x[, 1],
                pcab$x[, 2],
                xlab = paste("PC1:", pcabVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcabVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-batch"
        )
        
        plot(
                pcae$x[, 1],
                pcae$x[, 2],
                xlab = paste("PC1:", pcaeVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaeVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-error"
        )
        
        plot(
                pcac$x[, 1],
                pcac$x[, 2],
                xlab = paste("PC1:", pcacVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcacVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-corrected"
        )
}

#' svaplot is used to show the selected peaks as heatmap according to p-value and q-value
#' Input data could be the results from svacor function
#' Two methods could be used: pqvalues = "sv" means the surrogate variables analysis would be used to compute the correct data, p-value and q-value are also based on the corrected data; pqvalues = "annova" means the visulization of raw data according to p-value and q-value.
#' If there are seleced peaks to be show, svaplot will also return a list with selected data matrix
svaplot <- function(list,
                    pqvalues = "sv",
                    pt = 0.05,
                    qt = 0.05) {
        data <- list$data
        signal <- list$signal
        signal2 <- list$signal2
        batch <- list$batch
        error <- list$error
        error2 <- list$error2
        datacor <- list$dataCorrected
        pValues <- list$'p-values'
        qValues <- list$'q-values'
        pValuesSv <- list$'p-valuesCorrected'
        qValuesSv <- list$'q-valuesCorrected'
        
        icolors <-
                colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
        
        if (is.null(signal2)) {
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('No SV while p-values and q-values have results')
                        layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5), 9
                        ), 9, 9, byrow = TRUE))
                        par(mar = c(3, 5, 1, 1))
                        data <- data[pValues < pt & qValues < qt, ]
                        signal <-
                                signal[pValues < pt &
                                               qValues < qt, ]
                        error <-
                                error[pValues < pt & qValues < qt, ]
                        zlim <- range(c(data, signal, error))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt)
                        names(li) <- c('data', 'pqvalues')
                        return(li)
                }
                else{
                        message('No SV while p-values and q-values have no results')
                        layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4), 8
                        ), 8, 8, byrow = TRUE))
                        par(mar = c(3, 6, 2, 3))
                        zlim <- range(c(data, signal, error))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
        } else{
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('Have SVs while p-values and q-values have results')
                        layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5), 9
                        ), 9, 9, byrow = TRUE))
                        par(mar = c(3, 5, 1, 1))
                        data <- data[pValues < pt & qValues < qt, ]
                        signal <-
                                signal2[pValues < pt &
                                                qValues < qt, ]
                        error <-
                                error2[pValues < pt &
                                               qValues < qt, ]
                        zlim <- range(c(data, signal, error))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt)
                        names(li) <- c('data', 'pqvalues')
                        return(li)
                }
                else if (pqvalues == "anova") {
                        message('Have SVs while p-values and q-values have no results')
                        layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4), 8
                        ), 8, 8, byrow = TRUE))
                        par(mar = c(3, 6, 2, 3))
                        zlim <- range(c(data, signal2, error2))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(signal2),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(signal2) - 1
                                )),
                                labels = colnames(signal2),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(error2),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error2
                                ) - 1)),
                                labels = colnames(error2),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
                else if (pqvalues == "sv" &
                         sum(pValuesSv < pt &
                             qValuesSv < qt) != 0) {
                        message('SVs corrected while p-values and q-values have results')
                        layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6), 11
                        ), 11, 11, byrow = TRUE))
                        par(mar = c(3, 4, 2, 1))
                        data <-
                                data[pValuesSv < pt &
                                             qValuesSv < qt, ]
                        signal <- signal[pValuesSv < pt &
                                                 qValuesSv < qt, ]
                        batch <-
                                batch[pValuesSv < pt &
                                              qValuesSv < qt, ]
                        error <-
                                error[pValuesSv < pt &
                                              qValuesSv < qt, ]
                        datacor <-
                                datacor[pValuesSv < pt &
                                                qValuesSv < qt, ]
                        zlim <-
                                range(c(data, signal, batch, error, datacor))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(batch),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-batch',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        batch
                                ) - 1)),
                                labels = colnames(batch),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        batch
                                ) - 1)),
                                labels = rownames(batch),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(datacor),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-corrected',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(datacor) - 1
                                )),
                                labels = colnames(datacor),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (
                                        nrow(datacor) - 1
                                )),
                                labels = rownames(datacor),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(datacor,
                                     data,
                                     pValuesSv < pt &
                                             qValuesSv < qt)
                        names(li) <-
                                c('dataCorrected', 'data', 'pqvalues')
                        return(li)
                }
                else{
                        message('SVs corrected while p-values and q-values have no results')
                        layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6), 13
                        ), 13, 13, byrow = TRUE))
                        par(mar = c(3, 6, 2, 3))
                        zlim <-
                                range(c(signal, data, batch, error, datacor))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(batch),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-batch',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        batch
                                ) - 1)),
                                labels = colnames(batch),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 6, 2, 3))
                        image(
                                t(datacor),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-corrected',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(datacor) - 1
                                )),
                                labels = colnames(datacor),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
        }
}
#' Plot the comparison of PCA analysis before and after batch correction
svap <- function(name = 'pca.png',p1,p2,lv){
        png(name)
        par(mfrow = c(1,2))
        pca <- prcomp(t(p1))
        pcaVars=signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
        plot(pca$x[,1], 
             pca$x[,2], 
             xlab=paste("PC1:",pcaVars[1],"% of Variance Explained"),
             ylab=paste("PC2:",pcaVars[2],"% of Variance Explained"),
             pch=as.character(lv),
             cex=2,
             main = "PCA-Raw")
        pca <- prcomp(t(p2))
        pcaVars=signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
        plot(pca$x[,1], 
             pca$x[,2], 
             xlab=paste("PC1:",pcaVars[1],"% of Variance Explained"),
             ylab=paste("PC2:",pcaVars[2],"% of Variance Explained"),
             pch=as.character(lv),
             cex=2,
             main = "PCA")
        dev.off()
}
#' Plot the Quantitative visualization of Batch effects
plotb <- function(name = 'relativep.pdf', df, dfsv, dfanova,pos = 'topleft') {
        pdf(name, height = 8, width = 6)
        par(mfrow = c(2, 1))
        plot(
                df$PosteriorProbabilitiesSurrogate[dfanova$pqvalues] ~ df$PosteriorProbabilitiesMod[dfanova$pqvalues],
                xlab = 'Influences from experimental design',
                ylab = 'Influences from Batch effects',
                main = 'Influences of selected peaks(FDR control at 0.05)',
                cex = -log10(df$`p-valuesCorrected`[dfsv$pqvalues])-2,
                ylim = c(0, 1),
                xlim = c(0.6, 1.02)
        )
        points(
                df$PosteriorProbabilitiesSurrogate[dfsv$pqvalues] ~ df$PosteriorProbabilitiesMod[dfsv$pqvalues],
                pch = 1,
                cex = -log10(df$`p-values`[dfanova$pqvalues])-2,
                col = 'red'
        )
        legend(
                pos,
                c("raw", "corrected"),
                pch = c(1, 1),
                col = c('black', 'red')
        )
        plot(
                df$PosteriorProbabilitiesSurrogate ~ df$PosteriorProbabilitiesMod,
                xlab = 'Influences  from experimental design',
                ylab = 'Influences  from Batch effects',
                main = 'Influences  of all peaks',
                pch = 1,
                cex = 1.2,
                ylim = c(0, 1),
                xlim = c(0, 1.02)
        )
        dev.off()
}
# The following script is a demo to analysis metabolomics data in the main text
# Set data path
path <- "./pos/"
# Extrate data with optimized parameters for xcms
xset <- getdata(path,pmethod = 'uplcqtof')
# Assign the group information
lv <- as.factor(c(rep('C',10),rep('1',10),rep('2',10),rep('3',10)))
# Corrected the batch effects
df <- svacor(xset,lv)
# Plot the PCA analysis results in a linear model
svapca(df,center = F,scale = T,lv = as.character(lv))
# Plot the selected peaks as heatmap according to p-value and q-value
dfanova <- svaplot(df,pqvalues = 'anova')
dfsv <- svaplot(df,pqvalues = 'sv')
# Plot the PCA analysis before and after batch correction
svap(name = 'pca.png',dfanova$data,dfsv$dataCorrected,lv)
# plot the Quantitative visualization of Batch effects
plotb(name = 'relativep.pdf',df,dfsv,dfanova,pos = 'topleft')
