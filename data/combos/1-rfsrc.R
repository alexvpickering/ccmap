#R --max-ppsize 500000
library(randomForestSRC)
library(matrixStats)

#load data
data <- readRDS("~/Documents/Batcave/GEO/1-meta/data/MAMA/combos/combo_array.rds")$data

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
mv.obj1 <- rfsrc(f, data=data, ntree=800, nsplit=4, nodesize=1)
end.time <- Sys.time()

#how long?
time.taken <- end.time - start.time
time.taken

#error rates
err_tbl <- lapply(mv.obj1$classOutput, function(x) x$err.rate[, "all"])
err_tbl <- as.data.frame(err_tbl)

mean_err <- apply(err_tbl, 1, mean)
plot(1:length(mean_err), mean_err)











#compare OOB predictions to actual:
#----------------------------------
y <- data[, cb]

#get average OOB preds (averages across trees where sample is OOB)
oob1 <- lapply(mv.obj1$regrOutput, function(x) x$predicted.oob)
oob1 <- as.data.frame(oob1)

#compare sign of actual (y) to oob predictions (oob)
comp <- sign(oob2) == sign(y)

#per sample
comp_samples <- rowCounts(comp)
comp_samples <- comp_samples / ncol(y)



#combine OOB predictions from multiple RFs:
#------------------------------------------
oob2 <- lapply(mv.obj2$regrOutput, function(x) x$predicted.oob)
oob2 <- as.data.frame(oob2)

#get OOB counts for each sample
oob1_counts <- rowCounts(mv.obj1$inbag, value=0)
oob2_counts <- rowCounts(mv.obj2$inbag, value=0)

#combine using OOB counts to weight
w1 <- oob1_counts / (oob1_counts + oob2_counts)
w2 <- oob2_counts / (oob1_counts + oob2_counts)

oob12 <- (w1 * oob1) + (w2 * oob2)



#sample names
combo_data <- readRDS("/home/alex/Documents/Batcave/GEO/1-meta/combo_data.rds")
n_sels <- sapply(combo_data$selections, function(x) length(x))
sample_names <- rep(names(n_sels), n_sels)

#plot actual vs preds
for (i in 1:80){
    plot(as.numeric(oob1[i,]), 
         as.numeric(y[i,]), 
         main=sample_names[i], 
         xlab="pred", 
         ylab="actual")
}








