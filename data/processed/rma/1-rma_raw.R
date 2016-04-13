library(affy)
library(data.table)

setwd("~/Documents/Batcave/GEO/2-cmap/data/processed/rma")

#get cmap metadata
cmap_instances <- read.table("~/Documents/Batcave/GEO/2-cmap/data/raw/cmap_instances_02.csv", 
                             header=T, sep="\t", quote='', fill=T, stringsAsFactors=F)

cmap_instances <- data.table(cmap_instances)

#HT_HG-U133A_EA & HG-U133A: use below

array_instances <- cmap_instances[array3 =="HG-U133A",]

celfiles = c() # Each cel file has multiple controls
controls = list()  # So we keep them in a list

for (i in 1:length(array_instances$perturbation_scan_id)) {
  # This is the instance (cell exposed to drug) CEL file
  id = as.character(array_instances$perturbation_scan_id[i])
  id = gsub("'","",id)
  # This is the control cell (not exposed to drug)
  cid = as.character(array_instances$vehicle_scan_id4[i])
  
  #add path for instance file
  file = paste('~/Documents/Batcave/GEO/2-cmap/data/raw/',id,".CEL",sep="")
  celfiles = c(celfiles,file) 
  
  #add paths for control files (multiple controls)
  if (length(strsplit(cid, "[.]")[[1]]) > 2) {
    
    id = strsplit(id,"[.]")[[1]][1]
    cid = strsplit(cid,"[.]")
    cid = cid[[1]][-which(cid[[1]] %in% "")]
    
    tmp = c()
    for (c in 1:length(cid)){
      cinstance = paste(id,".",cid[c],sep="")
      file = paste('~/Documents/Batcave/GEO/2-cmap/data/raw/',cinstance,".CEL",sep="")
      tmp = c(tmp,file)    
    }
    controls = c(controls,list(tmp))
    
  #add paths for control files (single control)
  } else {
    tmp = c()
    for (c in 1:length(cid)){
      file = paste('~/Documents/Batcave/GEO/2-cmap/data/raw/',cid,".CEL",sep="")
      tmp = c(tmp,file)    
    }
    controls = c(controls,list(tmp)) 
    
  }
}

cel_paths <- unique(c(celfiles, unlist(controls)))
raw_data <- ReadAffy (filenames=cel_paths)
data <- affy::rma(raw_data)

saveRDS(data, "rma_HG-U133A.rds")


#HT_HG-U133A: 
#   - remove all "/home/alex/Documents/Batcave/GEO/2-cmap/data/raw/" & run up to cel_paths
#   - library(xps) #only worked in terminal R
#   - library(Biobase)
#   - scheme.hthgu133a.na35 <- root.scheme("/home/alex/Documents/Batcave/affy/schemes/na35/hthgu133a.root")
#   - celdir <- "/home/alex/Documents/Batcave/GEO/2-cmap/data/raw"
#   - cmap_hthgu133a <- import.data(scheme.hthgu133a.na35, "cmap_hthgu133a", celdir=celdir, celfiles=cel_paths, verbose=TRUE)
#   - cmap_hthgu133a_rma <- rma(cmap_hthgu133a, "cmap_hthgu133a_rma", verbose=FALSE)
#   - eset <- new("ExpressionSet", exprs = as.matrix(cmap_hthgu133a_rma))
#   - saveRDS(eset, "HT_HG-U133A.rds")
