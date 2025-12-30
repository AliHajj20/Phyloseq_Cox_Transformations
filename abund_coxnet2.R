abund_coxnet2 <- function(x_obj, time, status, covar = NULL, lambda = "lambda.1se", 
                         nvar = NULL, alpha = 0.9, nfolds = 10, showPlots = TRUE, 
                         coef_threshold = 0, transform = c("none", "logratio", "clr", "ilr", "alr", "rel")) {
  
  # Load necessary package
  library(compositions)
  
  # Transform argument to match one of the options
  transform <- match.arg(transform)
  
  # Impute zeros (optional, depending on data)
  x = data.frame(t(otu_table(x_obj)))
  x <- impute_zeros(x)
  
  # Set up survival object
  y = Surv(time, status)
  
  # Apply transformation to x
  if (transform == "logratio") {
    lrmatrix <- logratios_matrix(x)
    x_model <- lrmatrix[[1]]   # Log-ratio transformed data
    idlrX <- lrmatrix[[2]]     # Indices for log ratios
    nameslrX <- lrmatrix[[3]]  # Names of variables involved in log ratios
  } else if (transform == "clr") {
    x_model <- clr(x)  # Centered Log-Ratio transformation
  } else if (transform == "ilr") {
    x_model <- ilr(x)  # Isometric Log-Ratio transformation
    colnames(x_model) <- colnames(x)[-1]
  } else if (transform == "alr") {
    x_model <- alr(x)  # Additive Log-Ratio transformation
  } else if (transform == "rel") {
    rel_abund <- transform_sample_counts(x_obj, function(x) x / sum(x))
    x = data.frame(t(otu_table(rel_abund)))
    x_model <-as.matrix(x)   # Additive Log-Ratio transformation
  } else {
    x_model <- as.matrix(x)  # No transformation
  }
  
  # Perform the Cox regression with glmnet
  if (is.null(covar)) {
    cvfit <- glmnet::cv.glmnet(x_model, y, family = "cox", type.measure = "C", 
                               alpha = alpha, nfolds = nfolds, keep = TRUE)
  } else {
    df0 <- data.frame(as.matrix(y), covar)
    model0 <- coxph(Surv(time, status) ~ ., data = df0)
    x0 <- predict(model0)
    cvfit <- glmnet::cv.glmnet(x_model, y, family = "cox", type.measure = "C", 
                               nfolds = nfolds, alpha = alpha, keep = TRUE, offset = x0)
  }
  
  # Show plot of CV fit if requested
  if (showPlots == TRUE) {
    plot(cvfit)
  }
  
  # Select the best lambda based on number of variables or predefined lambda value
  if (!is.null(nvar)) {
    rowlasso <- max(which(cvfit$glmnet.fit$df <= nvar))
    lambda <- cvfit$glmnet.fit$lambda[rowlasso]
  }
  
  lambdavalue <- lambda
  if (is.character(lambda)) {
    if (lambda == "lambda.1se") 
      lambdavalue <- cvfit$lambda.1se
    if (lambda == "lambda.min") 
      lambdavalue <- cvfit$lambda.min
  }
  
  idrow <- max(which(cvfit$glmnet.fit$lambda >= lambdavalue))
  coeflr <- as.vector(coef(cvfit, s = lambda))
  
  # Calculate log-contrast for logratio transformation
  if (transform == "logratio") {
    coeflogcontrast <- rep(0, ncol(x))
    for (i in seq_along(coeflr)) {
      coeflogcontrast[idlrX[i, 1]] <- coeflogcontrast[idlrX[i, 1]] + coeflr[i]
      coeflogcontrast[idlrX[i, 2]] <- coeflogcontrast[idlrX[i, 2]] - coeflr[i]
    }
    varlogcontrast <- which(abs(coeflogcontrast) > coef_threshold)
    coeflogcontrast <- coeflogcontrast[varlogcontrast]
    names.select <- colnames(x)[varlogcontrast]
    logcontrast <- as.matrix(log(x)[, varlogcontrast]) %*% coeflogcontrast
  } else {
    varlogcontrast <- which(abs(coeflr) > coef_threshold)
    coeflogcontrast <- coeflr[varlogcontrast]
    names.select <- colnames(x_model)[varlogcontrast]
    logcontrast <- as.matrix(x_model[, varlogcontrast]) %*% coeflogcontrast
  }
  
  # Calculate predictions
  if (is.null(covar)) {
    predictions <- as.numeric(predict(cvfit, x_model, s = lambdavalue))
  } else {
    predictions <- as.numeric(predict(cvfit, x_model, s = lambdavalue, newoffset = x0))
  }
  
  # Normalize coefficients for interpretation
  coeflogcontrast <- 2 * coeflogcontrast / sum(abs(coeflogcontrast))
  
  # Compute the C-index for model performance
  if (length(varlogcontrast) == 0) {
    Cindex_signature <- 0.5
  } else {
    Cindex_signature <- glmnet::Cindex(pred = predictions, y)
  }
  
  # Get cross-validation statistics
  mcvCindex <- cvfit$cvm[idrow]
  sdcvCindex <- cvfit$cvsd[idrow]
  
  # Plot the risk score and signature plots
  plot1 <- NULL
  plot2 <- NULL
  if (length(varlogcontrast) > 0) {
    plot1 <- plot_riskscore(predictions, x, time, status, showPlots = showPlots)
    plot2 <- plot_signature(names.select, coeflogcontrast, showPlots = showPlots)
  } else {
    print("No variables are selected. The risk score plot and the signature plot are not displayed.")
  }
  
  # Return results
  results <- list(taxa.num = varlogcontrast, 
                  taxa.name = names.select, 
                  log_contrast_coefficients = coeflogcontrast, 
                  risk.score = predictions, 
                  apparent_Cindex = Cindex_signature, 
                  mean_cv_Cindex = mcvCindex, 
                  sd_cv_Cindex = sdcvCindex, 
                  risk_score_plot = plot1, 
                  signature_plot = plot2)
  
  return(results)
}
