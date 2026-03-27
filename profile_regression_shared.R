library(tidyverse)
library(sf)
library(sftime)
library(units)
library(PReMiuM)

rm(list = ls())

relevel.factor = function(x, others) {
  x = as.character(x)
  others = as.character(others)
  lev = unique(x)
  x[x %in% others] = 'Other'
  return(factor(x, levels = c(setdiff(lev, others) ,'Other')))
}

Kpc = 10 # number of principal components

load("sample_data.rda")


yModel = "Bernoulli"
f1pr = event ~ -1 + 
  slope_mean + east_mean + north_mean +
  profile_mean + planform_mean + TPI_mean + 
  VRM_mean +
  cos_month + sin_month #+ lith_code
f1pr = update(f1pr, 
              formula(paste0(".~.+", paste("PC", 1:Kpc, collapse = "+", sep = ""))))
mf.pr = rbind(
  model.frame(f1pr, presences), model.frame(f1pr, absences)
)
inp = list()
inp$profileNames = names(mf.pr)
inp$profileNames = setdiff(inp$profileNames, c("num_events", "event"))
inp$discreteCovs = names(mf.pr)[sapply(mf.pr, is.factor)]
inp$continuousCovs = setdiff(inp$profileNames, inp$discreteCovs)
if (length(inp$discreteCovs)==0) {
  inp$xModel = 'Normal'
} else if (length(inp$discreteCovs)>0) {
  inp$xModel = 'Mixed'
} else inp$xModel = 'Discrete'
inp$outcome = names(mf.pr)[1]


f1fx = event ~ -1 + lith_code 
mf.fx = rbind(
  model.frame(f1fx, presences), 
  model.frame(f1fx, absences)
)
mf.fx$lith_code = factor(mf.fx$lith_code)
table(mf.fx$lith_code, mf.fx$event)
mf.fx$lith_code = relevel.factor(mf.fx$lith_code, c('B', 'CM', 'E', 'Gd', 'M', 'Li', 'Mw', 'Pr'))
levels(mf.fx$lith_code)
mf.fx = data.frame(model.matrix(f1fx, mf.fx)[,-1])

inp$fixedEffectsNames = colnames(mf.fx)
inp$fixedEffectsNames = setdiff(inp$fixedEffectsNames, c("num_events", "event"))

mf.pr = mf.pr %>% bind_cols(dplyr::select(mf.fx, inp$fixedEffectsNames)) %>% 
  mutate(across(where(is.factor), as.double))


nburn = 10000
nrep = 15000
thin = 1

mod1 = profRegr(
  covNames = inp$profileNames, xModel = inp$xModel,
  discreteCovs = inp$discreteCovs, 
  continuousCovs = inp$continuousCovs,
  fixedEffectsNames = inp$fixedEffectsNames,
  outcome = inp$outcome, outcomeT = NA,
  output='mod5', 
  data = mf.pr,
  nSweeps=nrep, nBurn=nburn, nProgress=100, nFilter=thin,
  yModel = yModel, run = TRUE
)
dissimObj <- calcDissimilarityMatrix(mod1)

save(f1pr, f1fx, mf.pr, mf.fx, mod1, dissimObj, 
     file=paste0(mod1$fileStem ,"_profile_regression.RData"))


clusObj <- calcOptimalClustering(dissimObj, maxNClusters = 20)

riskProfileObj <- calcAvgRiskAndProfile(clusObj, includeFixedEffects=T)

clusterOrderObj <- plotRiskProfile(
  riskProfileObj, 
  paste0(mod1$fileStem, '_summary.png')
)
save(clusObj, riskProfileObj, clusterOrderObj, file=paste0(mod1$fileStem, "_risk_profile.rda"))


