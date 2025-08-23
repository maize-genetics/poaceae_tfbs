# Runs phylogenetic permulations for mixed models with asreml
# Charlie Hale, 2025.07.24

# input_df should have named columns corresponding to the assemblyID, orthogroup, response variable, predictor variable, and covariates
permulate <- function(input_df,
                      response_var_name,
                      predictor_var_name,
                      sp.tree,
                      Kphylo,
                      n,
                      threads = 1) {
  # Prep data
  predictor_names <- colnames(input_df)[4:ncol(input_df)]
  nPred <- length(predictor_names)
  nc <- 5 + 3*nPred
  results_mat <- matrix(NA, nrow = n+1, ncol = nc,
                        dimnames = list(NULL,
                                        c("OG","responseVar","focal_predVar","model","phylo_explained",
                                          paste0("waldStat_", predictor_names),
                                          paste0("coef_",      predictor_names),
                                          paste0("pval_log10_",predictor_names))))
  results_mat[,"OG"          ] <- input_df$OG[1]
  results_mat[,"responseVar" ] <- response_var_name
  results_mat[,"focal_predVar"] <- predictor_var_name
  results_mat[,"model"       ] <- c("empirical", sprintf("perm_%06d", 1:n))
  
  # Empirical fit
  results_mat[1,5:nc] <- run_model_asreml(input_df, Kphylo)
  
  # Simulate n times
  named_pred <- setNames(input_df[[predictor_var_name]], input_df$assemblyID)
  rate      <- geiger::ratematrix(sp.tree, named_pred)
  simulated <- geiger::sim.char(sp.tree, par = rate, nsim = n, model = "BM", root = 0)
  
  # Permulations (parallel if threads>1)
  sims <- seq_len(n)
  compute_one <- function(sim) {
    v <- simulated[,,sim]
    perm = sort(named_pred)[rank(v)]
    df2 <- input_df
    df2[[predictor_var_name]] <- perm
    run_model_asreml(df2, Kphylo)
  }
  
  if (threads > 1) {
    library(parallel)
    perm_rows <- mclapply(sims, compute_one,
                          mc.cores = threads,
                          mc.preschedule = FALSE)
  } else {
    perm_rows <- lapply(sims, compute_one)
  }
  results_mat[2:(n+1), 5:nc] <- do.call(rbind, perm_rows)
  
  results_mat
}

# Function to run phylogenetic association model using asreml
# Input_df should have the first three columns as assemblyID, OG, and response variable, then the other columns as the predictors
run_model_asreml <- function(input_df, Kphylo, modelFamily = "gaussian") {
  # fit model within tryCatch
  tryCatch({
    OG <- input_df[1,2] # Assuming second column is the orthogroup ID
    resp_var_name <- colnames(input_df)[3] # Third column name is the response variable
    pred_names <- colnames(input_df)[4:ncol(input_df)] # Predictors start from fourth column
    n_traits <- length(pred_names) # Number of predictors
    form  <- as.formula(paste(resp_var_name, "~", paste(pred_names, collapse = " + ")))
    if(modelFamily == "gaussian") {
      model <- asreml(
      fixed  = form,
      random = ~ vm(assemblyID, Kphylo),
      ai.sing = TRUE,
      data   = input_df
    )
    }
    else if(modelFamily == "binomial") {
      model <- asreml(
        fixed  = form,
        random = ~ vm(assemblyID, Kphylo),
        family = asreml::asr_binomial(),
        ai.sing = TRUE,
        data   = input_df
      )
    }
    
    summ  <- summary(model)
    wald  <- wald.asreml(model)
    phylo_var <- summ$varcomp$component[1]
    resid_var <- summ$varcomp$component[2]
    phylo_explained <- phylo_var / (phylo_var + resid_var)
    coeffs <- model$coefficients$fixed[2:(1 + n_traits)]
    waldStats <- wald[(2:(1 + n_traits)), 2]
    pvals_log10 <- -(pchisq(wald[2:(1 + n_traits), 3], df = 1, log.p = TRUE, lower.tail = FALSE) / log(10))
        
    # Ensure that the number of coefficients and p-values match the number of traits
    waldStats <- if (length(waldStats) == n_traits) waldStats else rep(NA_real_, n_traits)
    coeffs     <- if (length(coeffs) == n_traits) coeffs else rep(NA_real_, n_traits)
    pvals_log10 <- if (length(pvals_log10) == n_traits) pvals_log10 else rep(NA_real_, n_traits)

    c(phylo_explained, waldStats, coeffs, pvals_log10)
  }, error = function(e) {
    n_traits <- length(colnames(input_df)[4:ncol(input_df)])
    c(
      NA_real_,
      rep(NA_real_, n_traits),
      rep(NA_real_, n_traits),
      rep(NA_real_, n_traits)
    )
  })
}

