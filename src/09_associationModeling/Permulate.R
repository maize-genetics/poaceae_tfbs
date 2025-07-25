# Runs phylogenetic permulations for motif x orthogroup models
# Charlie Hale, 2025.07.24

# input_df should have named columns corresponding to the assemblyID, orthogroup, response variable, predictor variable, and covariates
permulate <- function(input_df, response_var_name, predictor_var_name, sp.tree, phyloK, n) {
    ## 1. Initialize matrix to store results
    nPredictors <- ncol(input_df) - 3 # Assuming first three columns are assemblyID, OG, and response variable
    results_mat <- matrix(NA, nrow = n + 1, ncol = 5 + 2 * nPredictors)
    colnames(results_mat) <- c("OG", "responseVar", "focal_predVar", "model", "phylo_explained", 
                               sprintf("coef_ePC%d", 1:nPredictors), 
                               sprintf("pval_log10_ePC%d", 1:nPredictors))
    results_mat[,"OG"] <- rep(input_df$OG[1], times = n + 1)
    results_mat[,"responseVar"] <- response_var_name
    results_mat[,"focal_predVar"] <- predictor_var_name
    results_mat[,"model"] <- c("empirical", sprintf("perm_%06d", 1:n))

    # Extract predictor and response variables
    named_pred <- input_df[,predictor_var_name]
    names(named_pred) <- input_df$assemblyID 
    named_resp <- input_df[,response_var_name]
    names(named_resp) <- input_df$assemblyID

    ## 1. Run model with empirical data
    results_mat[1,5:ncol(results_mat)] <- run_model_asreml(input_df, phyloK)

    ## 2. Running permulated models
    # Calculate evolutionary rate from the focal predictor variable  using geiger::rate.matrix
    rate <- geiger::ratematrix(sp.tree, named_pred)
    # Simulate trait along tree n times using geiger::sim.char with a Brownian motion model
    simulated <- geiger::sim.char(sp.tree, par = rate, nsim = n, model = "BM", root = 0)
    names(simulated) <- names(named_pred)
    
    # Loop through each simulated set and run model
    for(sim in 1:n) {
        # Extract the simulated data for this iteration
        sim_data <- simulated[, , sim]
        # Re-order the empirical data by the simulated data
        permulated_for_model = named_pred[order(sim_data)]
        # Run association model with permulated data
        perm_input_df <- input_df
        perm_input_df[, predictor_var_name] <- permulated_for_model # Sub in permulated values for focal predictor variable
        results_mat[sim + 1, 5:ncol(results_mat)] <- run_model_asreml(perm_input_df, phyloK)
    }

    # Return a matrix with results from empirical and permulated models
    return(results_mat)
}

# Function to run phylogenetic association model using asreml
# Input_df should have the first three columns as assemblyID, OG, and response variable, then the other columns as the predictors
run_model_asreml <- function(input_df, phyloK) {
  # fit model within tryCatch
  tryCatch({
    OG <- input_df[1,2] # Assuming second column is the orthogroup ID
    resp_var_name <- colnames(input_df)[3] # Third column name is the response variable
    pred_names <- colnames(input_df)[4:ncol(input_df)] # Predictors start from fourth column
    n_traits <- length(pred_names) # Number of predictors
    form  <- as.formula(paste(resp_var_name, "~", paste(pred_names, collapse = " + ")))
    model <- asreml(
      fixed  = form,
      random = ~ vm(assemblyID, phyloK),
      ai.sing = TRUE,
      data   = input_df
    )
    summ  <- summary(model)
    wald  <- wald.asreml(model)
    phylo_var <- summ$varcomp$component[1]
    resid_var <- summ$varcomp$component[2]
    phylo_explained <- phylo_var / (phylo_var + resid_var)
    coeffs <- model$coefficients$fixed[2:(1 + n_traits)]
    pvals_log10 <- -(pchisq(wald[2:(1 + n_traits), 3], df = 1, log.p = TRUE, lower.tail = FALSE) / log(10))

    coeffs     <- if (length(coeffs) == n_traits) coeffs else rep(NA_real_, n_traits)
    pvals_log10 <- if (length(pvals_log10) == n_traits) pvals_log10 else rep(NA_real_, n_traits)

    c(phylo_explained, coeffs, pvals_log10)
  }, error = function(e) {
    n_traits <- length(colnames(input_df)[4:ncol(input_df)])
    c(
      NA_real_,
      rep(NA_real_, n_traits),
      rep(NA_real_, n_traits)
    )
  })
}

