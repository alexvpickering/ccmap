library(xgboost)

setwd("~/Documents/Batcave/GEO/1-meta/data/MAMA/combos")


#--------
# SETUP
#--------

#load data
combo_probe <- readRDS("combo_probe.rds")$data


#get label (1 if up)
y <- ifelse(combo_probe$combo_dprime 0, 1, 0)

#remove combo info from X
cb_cols <- grepl("combo", colnames(combo_probe))
X <- combo_probe[, !cb_cols]

#get out-of-bag (oob) preds from rfsrc
model_rfsrc <- readRDS("1-model_rfsrc.rds")

oob <- lapply(model_rfsrc$classOutput, function(x) x$predicted.oob[, "1"])
oob <- t(as.data.frame(oob))  #n_probes (22215) * n_samples(80)
X[, "rfsrc_oob"] <- unlist(as.data.frame(oob))

#-------
# MODEL
#-------

#grid search
grid <- expand.grid(subsample=1, 
                     colsample_bytree=1,
                     max.depth=10)
ntrees <- 20000

#Build a xgb.DMatrix object
dtrain <- xgb.DMatrix(data=as.matrix(X), label=y)
hyps <- c()

for (i in nrow(grid)) {
 
   #Extract Parameters to test
   sub <- grid[["subsample"]][i]
   col <- grid[["colsample_bytree"]][i]
   dep <- grid[["max.depth"]][i]
 
   history <- xgb.cv(data=dtrain, nrounds=ntrees, nfold=5, 
                   objective="binary:logistic", eta=0.1,                               
                   subsample=sub, colsample_bytree=col, max.depth=dep,
                   early.stop.round=200, maximize=F)
 
   history <- as.data.frame(history)
 
   #Save error of the best iteration
   error <- min(history$test.error.mean)
   names(error) <- paste(sub, col, dep, sep="-")
 
   cat ("\nerror:", error, "sub-col-dep:", sub, col, dep, "\n\n")
   hyps <- c(hyps, error)
 }
