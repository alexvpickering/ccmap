#R --max-ppsize 500000
library(randomForestSRC)
library(matrixStats)

#load data
combo_array <- readRDS("~/Documents/Batcave/GEO/1-meta/data/MAMA/combos/combo_array.rds")
data <- combo_array$data

#-------
# SETUP
#-------

#remove combo columns that are not dprime values
#filt: T -> have combo_^(dprime) [remove]; F -> not combo_^(dprime) [keep]
filt <- grepl("combo_[^(dprime)]", colnames(data))
data <- data[, !filt]

#make formula
cols <- colnames(data)
cb <- grepl("combo", cols)

#for multivariate, formula is: cbind(y1, y2, ..., yi) ~ .
y <- paste(cols[cb], collapse=", ")
y <- paste("cbind(", y, ")", sep="")

f <- paste(y, ".", sep=" ~ ")
f <- as.formula(f)

#make combo columns binary factors (1 if up, otherwise 0)
data[, cb] <- ifelse(data[, cb] > 0, 1, 0)
data[, cb] <- lapply(data[, cb], factor)


#----------
# ANALYSIS
#----------

#8 trees / 11.35 mins * (19 * 60) mins ~= 800 trees
start.time <- Sys.time()
mv.obj <- rfsrc(f, data=data, ntree=800, nsplit=4, nodesize=1)
end.time <- Sys.time()

#how long?
time.taken <- end.time - start.time
time.taken

#error rates
err_tbl <- lapply(mv.obj$classOutput, function(x) x$err.rate[, "all"])
err_tbl <- as.data.frame(err_tbl)

mean_err <- apply(err_tbl, 1, mean)
plot(1:length(mean_err), mean_err)


model <- list(model=mv.obj, order=combo_array$order)
saveRDS(model, 'mod1.rds')