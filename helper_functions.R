################################################################################
#
# helper_functions.r
# Author: Amrit Singh
#
################################################################################

# load libraries
library(tidyverse);

################################################################################
#
# improve ggplot axes aesthetics
# Description: improve axes aesthetics
# Inputs:      sizeStripFont: size of facet label,
#              xAngle: angle of x-axis label,
#              hjust: horizontal alignement of x-axis label,
#              vjust: vertical alignement of x-axis label,
#              xSize: size of x-axis labels,
#              ySize: size of y-axis labels,
#              xAxisSize: size of x-axis heading,
#              yAxisSize: size of y-axis heading
# Outputs:     layer to add to a ggplot
#
################################################################################
customTheme = function (sizeStripFont = 10,
  xAngle = 0, hjust = 0.5, vjust = 0.5, xSize = 10, ySize = 10,
  xAxisSize = 10, yAxisSize = 10) {
  theme(strip.background = element_rect(colour = "black", fill = "white",
    size = 1), strip.text.x = element_text(size = sizeStripFont),
    strip.text.y = element_text(size = sizeStripFont), axis.text.x = element_text(angle = xAngle,
      hjust = hjust, vjust = vjust, size = xSize, color = "black"),
    axis.text.y = element_text(size = ySize, color = "black"),
    axis.title.x = element_text(size = xAxisSize, color = "black"),
    axis.title.y = element_text(size = yAxisSize, color = "black"),
    panel.background = element_rect(fill = "white", color = "black"))
}

################################################################################
#
# Heatmap of p-values
# Description: plots the p-values that depict
#              the association between principal components and
#              demographic variables
# Inputs:      pcs: principal components, matrix of n samples x k components
#              demo: demographic variables, n samples x p variables
# Outputs:     interactive ggplot heatmap
#
################################################################################
pcaPvalueHeatmap = function(pcs, demo, method="pearson"){
  result <- list()
  pvalHeatmap <- matrix(0, nc = ncol(demo), nr = ncol(pcs))
  rownames(pvalHeatmap) <- colnames(pcs)
  colnames(pvalHeatmap) <- colnames(demo)
  ## compute p-values beteween PCs and demographics
  for(i in 1:ncol(pcs)){
    for(j in 1:ncol(demo)){
      pvalHeatmap[i,j] <- cor.test(pcs[,i], demo[,j],method=method)$p.value
    }
  }
  result$pvalHeatmap <- pvalHeatmap


  pvalHeatmap[pvalHeatmap < 0.01] <- 0.01
  pvalHeatmap[pvalHeatmap > 0.1] <- 1
  pvalHeatmap[pvalHeatmap > 0.01 & pvalHeatmap < 0.05] <- 0.05
  pvalHeatmap[pvalHeatmap > 0.05 & pvalHeatmap < 0.1] <- 0.1
  pvalHeatmap[pvalHeatmap == "0.01"] <- "p < 0.01"
  pvalHeatmap[pvalHeatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalHeatmap[pvalHeatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalHeatmap[pvalHeatmap == "1"] <- "p > 0.10"
  pvals <- pvalHeatmap %>% as.data.frame %>% mutate(Variable = rownames(.)) %>%
    gather(Threshold, Value, -Variable) %>% mutate(Threshold = factor(Threshold,
      levels = unique(Threshold))) %>%
    mutate(Value = factor(Value, levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10")))
  result$pvals <- pvals
  p <- ggplot(pvals, aes(Threshold, Variable)) +
    geom_tile(aes(fill = Value), colour = "white") + scale_fill_manual(values = rev(brewer.pal(n = 8,
      name = "Blues")[c(2, 4, 6, 8)])) + customTheme(sizeStripFont = 10,
        xAngle = -40, hjust = 0, vjust = 1, xSize = 10, ySize = 10,
        xAxisSize = 10, yAxisSize = 10) + xlab("") + ylab("")
  result$p <- p
  result
}

################################################################################
#
# Heatmap of correlation between pcs and covariates
# Description: correlations that depict
#              the association between principal components and
#              demographic variables
# Inputs:      pcs: principal components, matrix of n samples x k components
#              demo: demographic variables, n samples x p variables
# Outputs:     correlation heatmap
#
################################################################################
pcaCorHeatmap = function(pcs, demo, method="pearson"){
  corHeatmap <- matrix(0, nc = ncol(demo), nr = ncol(pcs))
  rownames(corHeatmap) <- colnames(pcs)
  colnames(corHeatmap) <- colnames(demo)
  ## compute p-values beteween PCs and demographics
  for(i in 1:ncol(pcs)){
    for(j in 1:ncol(demo)){
      corHeatmap[i,j] <- cor(pcs[,i], demo[,j], method=method)
    }
  }

  corHeatmap
}

################################################################################
#
# modified fdrCsSAM() from the csSAM R-library
# Description:
# Inputs:
# Outputs:
#
################################################################################
fdrCsSAM2 = function (G, cc, y, n, numcell, numgene, rhat, nperms, alternative = "two.sided",
  standardize = TRUE, medianCenter = TRUE, logRm = FALSE, logBase = 2,
  nonNeg = FALSE)
{
  numgene = ncol(G)
  rhatperm <- array(dim = c(nperms, numcell, numgene))
  perm = list()
  for (i in 1:nperms) {
    o = sample(1:length(y))
    ystar = y[o]
    for (curset in levels(y)) {
      perm[[curset]] = csfit(cc[ystar == curset, ], G[ystar ==
          curset, ], logRm, logBase)
    }
    rhatperm[i, , ] = csSAM(perm[[1]]$ghat, perm[[1]]$se,
      n[1], perm[[2]]$ghat, perm[[2]]$se, n[2], standardize,
      medianCenter, nonNeg)
  }
  cutp.g = matrix(NA, nrow = numcell, ncol = 100)
  numcut = ncol(cutp.g)
  fdr.g = ncall.g = nperm.g <- array(dim = c(numcell, numcut))
  for (j in 1:numcell) cutp.g[j, ] = seq(0, max(abs(rhat[j,
    ])), length = 100)
  for (i in 1:numcut) {
    for (curcell in 1:numcell) {
      if (alternative == "two.sided") {
        fdr.g[curcell, i] = sum(abs(rhatperm[, curcell,
          ]) > cutp.g[curcell, i])/nperms/sum(abs(rhat[curcell,
            ]) > cutp.g[curcell, i])
        ncall.g[curcell, i] = sum(abs(rhat[curcell, ]) >
            cutp.g[curcell, i])
      }
      if (alternative == "greater") {
        fdr.g[curcell, i] = sum(rhatperm[, curcell, ] >
            cutp.g[curcell, i])/nperms/sum(rhat[curcell,
              ] > cutp.g[curcell, i])
        ncall.g[curcell, i] = sum(rhat[curcell, ] > cutp.g[curcell,
          i])
      }
      if (alternative == "less") {
        fdr.g[curcell, i] = sum(rhatperm[, curcell, ] <
            -cutp.g[curcell, i])/nperms/sum(rhat[curcell,
              ] < -cutp.g[curcell, i])
        ncall.g[curcell, i] = sum(rhat[curcell, ] < -cutp.g[curcell,
          i])
      }
    }
  }
  fdr.g = pmin(fdr.g, 1)
  for (j in 1:numcell) {
    fdr.g[j, ] = make.monotone(fdr.g[j, ])
  }
  return(list(fdr.g = fdr.g, rhatperm = rhatperm, cutp.g = cutp.g,
    rhat = rhat, ncall.g = ncall.g, alternative))
}



runSAM2 = function (G, y, s0.sam = NULL, stand.r = TRUE)
{
  if (stand.r == TRUE) {
    G = scale(t(G), center = apply(G, 1, median), scale = F)
  }
  if (is.null(s0.sam)) {
    s0.sam = quantile(ttest.func2(G, y)$sd, 0.5, na.rm = TRUE)
  }
  tt.sam = ttest.func2(G, y, s0 = s0.sam)$tt
  return(tt.sam)
}


ttest.func2 = function (x, y, s0 = 0, sd = NULL) {
  n1 <- sum(y == levels(y)[1])
  n2 <- sum(y == levels(y)[2])
  p <- nrow(x)
  m1 <- rowMeans(x[, y == levels(y)[1], drop = F])
  m2 <- rowMeans(x[, y == levels(y)[2], drop = F])
  if (is.null(sd)) {
    sd <- sqrt(((n2 - 1) * varr(x[, y == levels(y)[2]], meanx = m2) +
        (n1 - 1) * varr(x[, y == levels(y)[1]], meanx = m1)) * (1/n1 +
            1/n2)/(n1 + n2 - 2))
  }
  numer <- m2 - m1
  dif.obs <- (numer)/(sd + s0)
  return(list(tt = dif.obs, numer = numer, sd = sd))
}


fdrSAM2 = function (G, y, nperms, tt.sam, alternative = "two.sided")
{
  numgene = ncol(G)
  perm = list()
  ttstar.sam <- array(dim = c(nperms, numgene))
  for (i in 1:nperms) {
    o = sample(1:length(y))
    ystar = y[o]
    ttstar.sam[i, ] = runSAM2(G, ystar)
  }
  cutp.sam = seq(0, max(abs(tt.sam)), length = 100)
  sigGene.sam <- vector(mode = "numeric", numgene)
  fdr.sam = ncall.sam = rep(NA, 100)
  for (i in 1:100) {
    if (alternative == "two.sided") {
      fdr.sam[i] = (sum(abs(ttstar.sam) > cutp.sam[i])/nperms)/sum(abs(tt.sam) >
          cutp.sam[i])
      ncall.sam[i] = sum(abs(tt.sam) > cutp.sam[i])
      sigGene.sam[which(abs(tt.sam) > cutp.sam[i])] = fdr.sam[i]
    }
    if (alternative == "greater") {
      fdr.sam[i] = (sum(ttstar.sam > cutp.sam[i])/nperms)/sum(tt.sam >
          cutp.sam[i])
      ncall.sam[i] = sum(tt.sam > cutp.sam[i])
      sigGene.sam[which(tt.sam > cutp.sam[i])] = fdr.sam[i]
    }
    if (alternative == "less") {
      fdr.sam[i] = (sum(ttstar.sam < -cutp.sam[i])/nperms)/sum(tt.sam <
          -cutp.sam[i])
      ncall.sam[i] = sum(tt.sam < -cutp.sam[i])
      sigGene.sam[which(tt.sam < -cutp.sam[i])] = fdr.sam[i]
    }
  }
  fdr.sam = pmin(fdr.sam, 1)
  fdr.sam = make.monotone(fdr.sam)
  return(list(fdr.sam = fdr.sam, ncall.sam = ncall.sam, ttstar.sam = ttstar.sam,
    sigGene.sam = sigGene.sam))
}


plotCsSAM2 = function (csSAMdata, SAMdata, alternative, cellID, numcell)
{
  cellnames = cellID
  labs = paste("cell type ", as.character(1:numcell), sep = "")
  par(mfrow = c(3, 3))
  plot(SAMdata$ncall.sam, SAMdata$fdr.sam, xlab = "## called",
    ylab = "FDR", type = "l", log = "x", ylim = c(0, 1))
  title(paste("SAM", alternative))
  for (i in 1:numcell) {
    plot(csSAMdata$ncall.g[i, ], csSAMdata$fdr.g[i, ], xlab = "# called",
      ylab = "FDR", type = "l", log = "x", ylim = c(0,
        1))
    title(paste(cellnames[i], alternative))
  }
}
