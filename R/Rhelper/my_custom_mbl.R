my_custom_mbl <- function(Xr, Yr, Xu, Yu = NULL, 
                          k, 
                          two_step = FALSE,
                          cor_diss = 0.1, 
                          cor_k_range = c(10, 500), 
                          ws = 71,
                          spike = NULL,
                          method = local_fit_wapls(5, 25),
                          diss_method = "pca", 
                          pc_selection = list("opc", 20),
                          control = mbl_control(return_dissimilarity = FALSE, 
                                                validation_type = c("NNv"),
                                                tune_locally = FALSE,
                                                range_prediction_limits = FALSE,
                                                allow_parallel = TRUE), 
                          group = NULL,
                          center = TRUE, scale = FALSE, verbose = TRUE, 
                          documentation = character(), 
                          msc = TRUE, ...) {
  
  
  if (!is.null(group)) {
    group <- factor(group[!is.na(Yr)])
  }
  
  
  Xr <- Xr[!is.na(Yr),]
  Yr <- Yr[!is.na(Yr)]
  
  if (msc) {
    Xr_msc <- msc(Xr)
    Xu_msc <- msc(Xu, attr(Xr_msc, "Reference spectrum:"))
  } else {
    Xr_msc <- Xr
    Xu_msc <- Xu
  }
  
  
  if (verbose)
    cat(paste0(nrow(Xr), " observations with response values\n"))
  
  if (two_step) {    
    if (verbose)
      cat("Filtering reference observations based on the correlation dissimilarity...\n")
    
    
    knn <- search_neighbors(Xr_msc, 
                            Xu_msc, 
                            diss_method = "cor", 
                            ws = ws, 
                            k_diss = cor_diss, 
                            k_range = cor_k_range,
                            center = center, 
                            scale = scale)
    
    
    
    if (verbose)
      cat(paste0(length(knn$unique_neighbors), " observations selected\n"))
    
    # knn$unique_neighbors <- 1:nrow(Xr)
    
    if (!is.null(group)) {
      group <- factor(group[knn$unique_neighbors])
    }
    
    if (msc) {
      Xr_msc <- msc(Xr[knn$unique_neighbors,])
      Xu_msc <- msc(Xu, attr(Xr_msc,"Reference spectrum:"))
    } else {
      Xr_msc <- Xr[knn$unique_neighbors,]
      Xu_msc <- Xu
    }
    
    Yr <- Yr[knn$unique_neighbors]
  }
  # mbl (same preprocessing as general models)
  mbl_results <- mbl(
    Xr_msc,
    Yr,
    Xu_msc, 
    Yu = Yu, 
    k = k,
    spike = NULL,
    method = method,
    diss_method = diss_method, 
    diss_usage = "none",
    gh = TRUE, 
    ws = ws,
    pc_selection = pc_selection,
    control = control, 
    group = group,
    seed = 179,
    center = center, 
    scale = scale, 
    verbose = verbose, 
    documentation = documentation)
  
  
  if (two_step)
    mbl_results$observations_used <- knn$unique_neighbors
  
  mbl_results$rmse_index_optimal <- which.min(mbl_results$validation_results$nearest_neighbor_validation$rmse)
  mbl_results
}


plot_my_custom_mbl <- function(x, what = "rmse", ...) {
  plot_dots <- list(...)
  preds <- as.data.frame(x$validation_results$Yu_prediction_statistics)
  nnval <- as.data.frame(x$validation_results$nearest_neighbor_validation)
  
  tpl <- cbind(preds[, c("k", what)],
               nnval[, c(what), drop = F])
  
  colnames(tpl) <- c("k", "Prediction", "NN validation")
  matplot(x = tpl[, 1], 
          y = tpl[, -1], 
          type = "b", 
          pch = 16,
          xlab = "k", 
          ylab = what, 
          ylim = c(min(tpl[, -1]), 1.1 * max(tpl[, -1])), 
          col = c("red", "dodgerblue"),
          ...)
  grid(nx = NULL, ny = NULL, col = rgb(0.3, 0.3, 0.3, 
                                       0.1), lty = 1, lwd = 1, equilogs = TRUE)
  legend("topright", 
         legend = colnames(tpl[,-1, drop = FALSE]), 
         bg = NA, 
         pch = 16, 
         col =  c("red", "dodgerblue"), 
         box.lty = 0)
  
  
}

plot_val <- function(pred, ref, label_pred = NULL, label_ref = NULL) {
  # cor(pred, ref)^2
  plot(pred,  
       ref, 
       xlim = range(ref, pred),
       ylim = range(ref, pred),
       xlab = label_pred,
       ylab = label_ref)
  abline(0, 1, col = "red")
  grid()
}

