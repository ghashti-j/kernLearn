#######################################################
# The MSCV objective function for optimal bandwidth calculation as introduced in
# Ghashti, J. S., & Thompson, J. R. (2023). Kernel Metric Learning for Clustering Mixed-type Data. arXiv preprint arXiv:2306.01890.

# as called in the kdsum() function from kdsum.R
# may also be used to calculated maximum similarity cross-validation bandwidth values based on the kdsum.R function

# To be used prior to using the kdsum function in the kdsum.R file as the bw argument vector
# Note that any bandwidth selection technique may be used, but will not have the features outline in above paper
#######################################################

mscv_dkss <- function(df, nstart = NULL, ckernel = "c_gaussian",
                       ukernel = "u_aitken", okernel = "o_wangvanryzin", verbose = TRUE) {
  if (is.null(nstart)) {
    nstart <- ifelse(ncol(df) > 4, 3, ncol(df))
    message("No nstart value given, defaulting to ", nstart)
  }
  v_ck <- c("c_gaussian", "c_epanechnikov", "c_uniform", "c_triangle",
            "c_biweight", "c_triweight", "c_tricube", "c_cosine", 
            "c_logistic", "c_sigmoid", "c_silverman")
  v_uk <- c("u_aitken", "u.aitchisonaitken")
  v_ok <- c("o_wangvanryzin", "o_habbema", "o_aitken", "o_aitchisonaitken", "o_liracine")
  
  if (!(ckernel %in% v_ck)) {
    stop("Invalid ckernel specified. Choose one of: ", paste(v_ck, collapse = ", "))
  }
  if (!(ukernel %in% v_uk)) {
    stop("Invalid ukernel specified. Choose one of: ", paste(v_uk, collapse = ", "))
  }
  if (!(okernel %in% v_ok)) {
    stop("Invalid okernel specified. Choose one of: ", paste(v_ok, collapse = ", "))
  }
  
  # Get column names by type
  con_cols <- names(df)[sapply(df, is.numeric)]
  fac_cols <- names(df)[sapply(df, function(x) is.factor(x) & !is.ordered(x))]
  ord_cols <- names(df)[sapply(df, is.ordered)]
  
  # Reorder the dataframe
  df_ordered <- df[, c(con_cols, fac_cols, ord_cols)]
  
  # store indices
  con_ind <- length(con_cols) #index of continuous variables
  fac_ind <- con_ind + length(fac_cols) #index of factors
  ord_ind <- fac_ind + length(ord_cols) #index of ordinal variables
  
  df_ordered <- data.matrix(df_ordered)
  n <- nrow(df_ordered)
  N <- ncol(df_ordered)
  
  # MSCV optimization function
  mscv_opt <- function(lambda) {
    penalty <- sum(c(
      if (con_ind > 0) lambda[1:con_ind] <= 0 else numeric(),
      if (fac_ind > con_ind) sapply((con_ind + 1):fac_ind, function(i) {
        max_val <- if (ukernel == "u.aitchisonaitken") {
          (max(df_ordered[, i]) - 1) / max(df_ordered[, i])
        } else 1
        lambda[i] < 0 || lambda[i] > max_val
      }) else numeric(),
      if (ord_ind > fac_ind) lambda[(fac_ind + 1):ord_ind] < 0 | lambda[(fac_ind + 1):ord_ind] > 1 else numeric()
    ))
    if (penalty > 0) return(Inf)
    
    outK <- outL <- outell <- K <-  L <-  ell <- list()
    
    if(con_ind != 0){
      for (i in 1:con_ind) { outK[[i]] <- outer(df_ordered[, i], df_ordered[, i], "-") }
    }
    if(fac_ind > con_ind){
      for (i in (con_ind + 1):fac_ind) { outL[[i - con_ind]] <- outer(df_ordered[, i], df_ordered[, i], "-") }
      outL <- Filter(Negate(is.null), outL)
    }
    if(ord_ind > fac_ind){
      for (i in (fac_ind + 1):ord_ind) { outell[[i - fac_ind]] <- outer(df_ordered[, i], df_ordered[, i], "-") }
      outell <- Filter(Negate(is.null), outell)
    }
    
    # Kernel calculations for continuous variables
    if (con_ind != 0) {
      K <- lapply(1:con_ind, function(i) {
        out <- outK[[i]]
        lambda_val <- lambda[i]
        switch(ckernel,
               "c_gaussian" = (1/sqrt(2*pi)) * exp(-0.5 * (out / lambda_val)^2),
               "c_epanechnikov" = (3/4) * (1 - (out / lambda_val)^2) * (abs(out / lambda_val) <= 1),
               "c_uniform" = 0.5 * (abs(out / lambda_val) <= 1),
               "c_triangle" = (1 - abs(out / lambda_val)) * (abs(out / lambda_val) <= 1),
               "c_biweight" = (15/16) * (1 - (out / lambda_val)^2)^2 * (abs(out / lambda_val) <= 1),
               "c_triweight" = (35/32) * (1 - (out / lambda_val)^2)^3 * (abs(out / lambda_val) <= 1),
               "c_tricube" = (70/81) * (1 - abs(out / lambda_val)^3)^3 * (abs(out / lambda_val) <= 1),
               "c_cosine" = (pi/4) * cos((pi/2) * out / lambda_val) * (abs(out / lambda_val) <= 1),
               "c_logistic" = 1 / (exp(out / lambda_val) + 2 + exp(-out / lambda_val)),
               "c_sigmoid" = 2 / (pi * (exp(out / lambda_val) + exp(-out / lambda_val))),
               "c_silverman" = 0.5 * exp(-abs(out) / sqrt(2)) * sin(abs(out) / sqrt(2) + pi/4)
        )
      })
      D1 <- Reduce(`+`, K)
    } else {
      D1 <- 0
    }
    
    # Kernel calculations for unordered factors
    if (fac_ind > con_ind) {
      L <- lapply(1:length(outL), function(i) {
        out <- outL[[i]]
        lambda_val <- lambda[i + con_ind]
        switch(ukernel,
               "u_aitken" = ifelse(out == 0, 1, lambda_val),
               "u.aitchisonaitken" = ifelse(out == 0, 1 - lambda_val, lambda_val / ((max(unique(out)) - 1)))
        )
      })
      D2 <- Reduce(`+`, L)
    } else {
      D2 <- 0
    }
    
    # Kernel calculations for ordered factors
    if (ord_ind > fac_ind) {
      ell <- lapply(1:length(outell), function(i) {
        out <- outell[[i]]
        lambda_val <- lambda[i + fac_ind]
        switch(okernel,
               "o_habbema" = lambda_val^(abs(out)^2),
               "o_wangvanryzin" = ifelse(out == 0, 1 - lambda_val, 0.5 * (1 - lambda_val) * lambda_val^abs(out)),
               "o_aitken" = ifelse(out == 0, lambda_val, (1 - lambda_val) / (2^abs(out))),
               "o_aitchisonaitken" = choose(max(unique(out)), abs(out)) * lambda_val^abs(out) * (1 - lambda_val)^(max(unique(out)) - abs(out)),
               "o_liracine" = ifelse(out == 0, 1, lambda_val^abs(out))
        )
      })
      D3 <- Reduce(`+`, ell)
    } else {
      D3 <- 0
    }
    
    D <- D1 + D2 + D3
    diag(D) <- 0
    Fx <- mean(log(colSums(D)) - log(n-1))
    return(Fx)
  }
  max_val <- -Inf
  best_obj <- NULL
  max_fail <- 7
  fail_cnt <- 0
  for (i in 1:nstart) {
    if (verbose == TRUE) print(paste0("start ", i, " of ", nstart))
    params <- runif(N, 1e-16, 1)
    while (TRUE) {
      tryCatch({
        result <- optim(par = params, mscv_opt, control = list(fnscale = -1, maxit = 10000))
        break  
      }, error = function(e) {
        fail_cnt <<- fail_cnt + 1
        if (fail_cnt >= max_fail) {
          stop("Too many optimization failures, please try different kernel functions or transforming variables that may be on a different scale")
        }
        print("Optimization failed, trying again...")
      })
    }
    
    if (result$value > max_val) {
      max_val <- result$value
      best_obj <- result
    }
  }
  if(verbose == TRUE) ifelse(best_obj$convergence == 0, print("Objective converged"), print("Objective did not converge."))
  bw <- data.frame(x = round(best_obj$par, 6))
  rownames(bw) <- colnames(df[, c(con_cols, fac_cols, ord_cols)])
  return(list(bw = bw, fn_value = best_obj$value))
}
