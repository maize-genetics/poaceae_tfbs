# Helper function to load motif counts into a list of matrices in R
load_motif_data <- function(unenriched_clusters_path, cluster_meta_path, motifs_dir) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(parallel)
  
  # Load unenriched clusters
  unenriched_clusters <- read.table(unenriched_clusters_path, header = FALSE)
  
  # Load motif cluster metadata and filter
  cluster.meta <- read.table(cluster_meta_path) %>%
    dplyr::select(TF = V3, cluster = V2) %>%
    dplyr::filter(!TF %in% unenriched_clusters$V1) %>%
    arrange(TF)
  
  # List motif files
  files <- list.files(motifs_dir, full.names = TRUE)
  n_ortho <- sum(file.info(files)$isdir == FALSE)
  
  # Parallel processing
  # Parallel processing
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = c("files", "n_ortho", "cluster.meta"), envir = environment())
  
  mat.list <- parLapply(cl, 1:n_ortho, function(i) {
    library(data.table)
    library(dplyr)
    library(tidyr)
    
    motifCounts <- fread(files[[i]]) %>%
      dplyr::select(c(cluster = V2, ogID = V3, motifCount = V4, assemblyID = V5)) %>%
      dplyr::left_join(cluster.meta) %>%
      dplyr::select(-cluster) %>%
      dplyr::distinct(across(c(assemblyID, TF)), .keep_all = TRUE) %>%
      pivot_wider(names_from = TF, values_from = motifCount)
    
    # Handle missing TFs
    missing_TFs <- setdiff(cluster.meta$TF, colnames(motifCounts))
    motifCounts[missing_TFs] <- as.integer(NA)
    motifCounts <- motifCounts %>%
      dplyr::select(assemblyID, ogID, all_of(cluster.meta$TF)) %>%
      dplyr::mutate(across(where(is.numeric), ~replace_na(., 0)))
    
    mat <- as.matrix(motifCounts[, 3:ncol(motifCounts)])
    rownames(mat) <- motifCounts$assemblyID
    return(mat)
  })
  
  stopCluster(cl)
  names(mat.list) <- gsub(".*/([^/]+)\\.txt", "\\1", files) # Annotate with OG IDs
  
  return(mat.list)
}
