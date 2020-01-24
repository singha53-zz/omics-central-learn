################################################################################
#
# helper_functions.r
# Author: Amrit Singh
#
################################################################################

# load libraries
library(tidyverse);
library(amritr);

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
  pvalHeatmap <- matrix(0, nc = ncol(demo), nr = ncol(pcs))
  rownames(pvalHeatmap) <- colnames(pcs)
  colnames(pvalHeatmap) <- colnames(demo)
  ## compute p-values beteween PCs and demographics
  for(i in 1:ncol(pcs)){
    for(j in 1:ncol(demo)){
      pvalHeatmap[i,j] <- cor.test(pcs[,i], demo[,j],method=method)$p.value
    }
  }

  pvalHeatmap[pvalHeatmap < 0.01] <- 0.01
  pvalHeatmap[pvalHeatmap > 0.1] <- 1
  pvalHeatmap[pvalHeatmap > 0.01 & pvalHeatmap < 0.05] <- 0.05
  pvalHeatmap[pvalHeatmap > 0.05 & pvalHeatmap < 0.1] <- 0.1
  pvalHeatmap[pvalHeatmap == "0.01"] <- "p < 0.01"
  pvalHeatmap[pvalHeatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalHeatmap[pvalHeatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalHeatmap[pvalHeatmap == "1"] <- "p > 0.10"
  pvalHeatmap %>% as.data.frame %>% mutate(Variable = rownames(.)) %>%
    gather(Threshold, Value, -Variable) %>% mutate(Threshold = factor(Threshold,
      levels = unique(Threshold))) %>%
    mutate(Value = factor(Value, levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10"))) %>%
    ggplot(aes(Threshold, Variable)) +
    geom_tile(aes(fill = Value), colour = "white") + scale_fill_manual(values = rev(brewer.pal(n = 8,
      name = "Blues")[c(2, 4, 6, 8)])) + customTheme(sizeStripFont = 10,
        xAngle = 40, hjust = 1, vjust = 1, xSize = 10, ySize = 10,
        xAxisSize = 10, yAxisSize = 10) + xlab("") + ylab("")
}
