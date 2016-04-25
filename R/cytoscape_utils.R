

get_tan_sim <- function(drug_info, es=F) {
  #implements drug similarity measurement from DMAP (Tanimoto Coefficient)
  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4597058/#B11

  #get top up/dn drug genes
  drug_genes <- get_drug_genes(drug_info, drug_genes_n=500, es=es)

  #compute net similarity between each drug
  sim <- list()
  for (i in seq_along(drug_genes)) {

    #determine direct overlap
    query_genes <- drug_genes[[i]]
    overlap <- get_overlap(drug_genes, query_genes, tan_sim=T)

    #get cross overlap (drug acts opposite to query)
    cross_genes <- list(up=query_genes$dn, dn=query_genes$up)
    cross_overlap <- get_overlap(drug_genes, cross_genes, tan_sim=T)

    #compute net overlap
    overlap$table <- get_net_overlap(overlap, cross_overlap, tan_sim=T)

    #add overlap
    sim[[names(drug_genes)[i]]] <- overlap
  }

#get similarity matrix
sim_mat <- get_tan_sim_mat(sim)

return(sim_mat)
}

#---------------------

get_tan_sim_mat <- function(sim) {
  #converts sim from tan_sim into correlation matrix

    drug_names <- names(sim)
    sim_mat <- matrix(nrow=length(drug_names),
                      ncol=length(drug_names),
                      dimnames=list(drug_names, drug_names))

    for (dr in drug_names){
        drug_table <- sim[[dr]]$table
        sim_mat[row.names(drug_table), dr] <- drug_table$net
    }
    return(sim_mat)
}


#---------------------

get_gsea_sim <- function(drug_info) {
  #setup ExpressionSet for similarity measuring
  pheno_data <- data.frame(state = colnames(drug_info),
                           row.names = colnames(drug_info))

  pheno_data <- as(pheno_data, "AnnotatedDataFrame")
  drug_info_eset <- ExpressionSet(drug_info, phenoData = pheno_data)

  #distance/similarity measuring
  ds <- GeneExpressionSignature::ScoreGSEA(drug_info_eset, 250, "avg")
  sim_mat <- (1-ds)

  return(sim_mat)
}

#---------------------

make_cytofiles <- function(sim_mat) {

  #affinity propogation clustering
  cluster_result <- apcluster(sim_mat)

  #generate network/similarity file for cytoscape
  exemplars <- names(cluster_result@exemplars)
  clusters <- lapply(cluster_result@clusters, names)

  network <- data.frame(Source = character(0),
                        Type = character(0),
                        Target = character(0), stringsAsFactors = F)

  sim <- data.frame(Edge = character(0),
                    Similarity = numeric(0), stringsAsFactors = F)

  for (cluster in clusters) {
    #add value of closest exemplar
    exemplar <- intersect(cluster, exemplars)  #source
    closest <- names(which.max(sim_mat[setdiff(exemplars, exemplar), exemplar]))  #target
    value <- sim_mat[closest, exemplar]

    edge <- paste(exemplar, "(1)", closest)
    row <- nrow(network) + 1

    network[row, ] <- c(exemplar, "1", closest)
    sim[row, ] <- c(edge, value)

    for (drug in setdiff(cluster, exemplar)) {
      #add value of cluster drugs
      value <- sim_mat[drug, exemplar]

      edge <- paste(drug, "(1)", exemplar)
      row <- nrow(network) + 1

      network[row, ] <- c(drug, "1", exemplar)
      sim[row, ] <- c(edge, value)
    }
  }

  #save network/similarity files
  write.table(network, "drug_network.txt", quote=FALSE, row.names=FALSE, sep='\t')
  write.table(sim, "drug_sim.txt", quote=FALSE, row.names=FALSE, sep='\t')

  #generate drug clusters file
  drug_names <- colnames(sim_mat)
  drug_clusters <- data.frame(Drug=drug_names)

  for (i in seq_along(cluster_result)) {
    drug_clusters[cluster_result[[i]], "Cluster"] <- i
  }
  write.table(drug_clusters, "drug_clusters.txt", quote=FALSE, row.names=FALSE, sep='\t')
}

#-------------------
