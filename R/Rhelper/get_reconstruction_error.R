# pls reconstruction analyisis
# this function takes Xr and Yr to build a pls projection model and
# then it uses this model to project Xu. The Xu matrix is 
# back-transformed/reconstructed and compared to the original Xu to finally
# obtain the reconstruction errors
# Arguments:
# Xr afsis spectra
# Yr afsis response variable(s). More than one are accepted
# Xu congo spectra
# pc_selection how to choose the number of components. See ortho_projection 
#              function in package resemble
# Return:
# error The reconstruction error
# scaled_error The reconstruction error scaled based on the range of values in Xu
# npls The number of components used

get_reconstruction_error <- function(Xr, 
                                     Yr, 
                                     Xu, 
                                     pc_selection = list("opc", 20)) {
  
  pls_proj <- ortho_projection(Xr,
                               Yr = Yr,
                               method = "pls",
                               pc_selection = pc_selection,
                               scale = TRUE)

  congo_proj <- predict(pls_proj, Xu)
  
  # back-transform (reconstruct) to spectral space using the pls loadings
  congo_back_transformed <- congo_proj %*% (pls_proj$X_loadings)
  
  scaled_xu <- scale(Xu, TRUE, TRUE)
  
  # finish back-transform by back-scaling 
  congo_back_transformed <- sweep(congo_back_transformed,
                                  MARGIN = 2,
                                  FUN = "*",
                                  STATS = as.vector(pls_proj$scale))
  
  congo_back_transformed <- sweep(congo_back_transformed,
                                  MARGIN = 2,
                                  FUN = "+",
                                  STATS = colMeans(Xr))
  
  # compute reconstruction error
  rec_error <- mean((scaled_xu - scale(congo_back_transformed, TRUE, TRUE))^2)^0.5
  # plot(scaled_xu[1,], type = "l")
  # plot(scale(congo_back_transformed, T, T)[1,], type = "l")
  
  
  # scale error for comparisons with other results
  rec_error_scaled <- rec_error/diff(range(scaled_xu))
  
  
  eval_result <- list(error = rec_error,
                      scaled_error = rec_error_scaled,
                      npls = pls_proj$n_components)
  
  eval_result <- c(eval_result, 
                   as.list(pls_proj$variance$y_var[, pls_proj$n_components]))
  
  if ("opc" %in% pc_selection) {
    eval_result[["rmsd"]] <- min(pls_proj$opc_evaluation[,ncol(pls_proj$opc_evaluation)])  
  }
  eval_result
}
