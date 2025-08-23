# Calculates motif enrichment in empirical vs background regions.
# Charlie Hale, 2025.07.10
# File names should begin with the motif ID, e.g. "MA0001.1.fimo.txt"
calculateMotifEnrichment <- function(motif_id, emp_dir, bg_dir, method = "counts") {
    library(fitdistrplus)
    library(tibble)
    library(data.table)
    library(purrr)
    options(scipen = 20)

    # Load empirical motif data
    emp_file <- list.files(emp_dir, pattern = paste0("^", motif_id), full.names = TRUE)[1] # Empirical regions
    bg_files  <- list.files(bg_dir,  pattern = paste0("^", motif_id), full.names = TRUE) # Background regions
    if ( method == "scores" ) {
        # Enrichment based on motif scores
        emp <- sum(data.table::fread(emp_file)$score)
        bg <- lapply(bg_files, function(file) sum(data.table::fread(file)$score)) %>%
            unlist()
    } else if (method == "counts") {
        # Enrichment based on motif counts
        emp <- nrow(data.table::fread(emp_file))
        bg <- lapply(bg_files, function(file) nrow(data.table::fread(file))) %>%
            unlist()
    }
    
    # Score-based enrichment
    fit <- fitdistr(bg, "normal")
    fc <- emp / fit$estimate["mean"]
    pvalue <-   2 * pnorm(-abs(emp - fit$estimate["mean"]),
                mean = 0,
                sd = fit$estimate["sd"])
    # Return a vector with motif ID, empirical score, background mean, background sd, and two-tailed p-value
    c(motif_id, emp, fit$estimate["mean"], fit$estimate["sd"], pvalue, fc, method)
}