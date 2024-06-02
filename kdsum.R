###########################################
# Computes the kernel summation distance for the kdsum function presented in 
# Ghashti, J. S., & Thompson, J. R. (2023). Kernel Metric Learning for Clustering Mixed-type Data. arXiv preprint arXiv:2306.01890.

# function(cIND, uIND, oIND, cFUN = "c_gaussian", uFUN = "u_aitken", oFUN = "o_wangvanryzin", bw, df, stan = FALSE)
# see readME file for citations on kernels used, as well as formulas. Kernel function options can be seen in lines 27--44 below

# Assuming the dataframe as ordered as continuous variable first, then nominal, then ordinal variable --
# cIND: the index of the last continuous variable in the dataframe (enter 0 if no continuous variables)
# uIND: the index of the last nominal variable in the dataframe (enter 0 if no nominal variables)
# oIND: the index of the last ordina variable in the dataframe (enter 0 if no ordinal variables)
# cFUN: the kernel function to be used for the continuous variables (enter NA if no continuous variables); defaults to the Gaussian/RBF kernel
# uFUN: the kernel function to be used for the nominal variables (enter NA if no nominal variables); defaults to the Aitken Kernel
# oFUN: the kernel function to be used for the ordinal variables (enter NA if no ordinal variables); defaults to the Wang & van Ryzin kernel
# bw: a vector of bandwidths, likely obtained from the MSCV objective function (may use any other bandwidths). Must be same length as total number of variables in the dataframe.
# stan: TRUE if the distance matrix should be standardized between [0,1]; FALSE otherwise

## INPUT: a dataframe
## OUTPUT: a matrix of pairwise distances
#############################################

kdsum <- function(cIND, uIND, oIND, cFUN = "c_gaussian", uFUN = "u_aitken", oFUN = "o_wangvanryzin", bw, df, stan = FALSE) {
  df <- data.matrix(df)
  if(length(bw) != ncol(df)) {stop("Invalid length of bandwidth vector (bw).")}
  if(any(bw[1:cIND] == 0)) {stop("Continuous bandwidths must be > 0")}
  if(uIND > 0 && any(bw[(cIND + 1):uIND] > 1)) {warning("Nominal bandwidths should be between 0 and 1")} 
  if(oIND > 0 && any(bw[(cIND + uIND + 1):oIND] > 1)) {warning("Ordinal bandwidths should be between 0 and 1")} 
  distances <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  combinations <- combn(nrow(df), 2)
  kernel_functions <- list(
    c_gaussian = function(A, B, bw, C) {sum((2 * pi)^(-1/2) * exp(-((A - B) / bw)^2 / 2))}, #Continuous Gaussian kernel
    c_uniform = function(A, B, bw, C) {sum(sapply(1:length(A), function(i) ifelse(abs((A[i] - B[i]) / bw[i]) <= 1, 1/2, 0)))}, #Continuous Uniform kernel
    c_epanechnikov = function(A, B, bw, C) {sum(ifelse(abs((A - B) / bw) <= 1, (3/4) * (1 - ((A - B) / bw)^2), 0))}, #Continuous Epanechnikov kernel
    c_triangle = function(A, B, bw, C) {sum(ifelse(abs((A - B) / bw) <= 1, (1 - abs((A - B) / bw)), 0))}, #Continuous Triangle kernel
    c_biweight = function(A, B, bw, C) {sum(ifelse(abs((A - B) / bw) <= 1, (15/32) * (3 - ((A - B) / bw)^2)^2, 0))}, #Continuous biweight kernel
    c_triweight = function(A, B, bw, C) {sum(ifelse(abs((A - B) / bw) <= 1, (35/32) * (1 - ((A - B) / bw)^2)^3, 0))}, #Continuous triweight kernel
    c_tricube = function(A, B, bw, C) {sum(ifelse(abs((A - B) / bw) <= 1, (70/81) * (1 - abs((A - B) / bw)^3)^3, 0))}, #Continuous tricube kernel
    c_cosine = function(A, B, bw, C) {sum(ifelse(abs((A - B) / bw) <= 1, (pi/4) * cos((pi/2) * ((A - B) / bw)), 0))}, #Continuous cosine kernel
    c_logistic = function(A, B, bw, C) {sum(1 / (exp((A - B) / bw) + 2 + exp(-((A - B) / bw))))}, #Continuous logistic kernel
    c_sigmoid = function(A, B, bw, C) {sum(2 / (pi * (exp((A - B) / bw) + exp(-((A - B) / bw)))))}, #Continuous Sigmoid kernel
    u_aitchisonaitken = function(A, B, bw, C) {sum(ifelse(A == B, 1 - bw, bw / (length(unique(C[,i])) - 1)))}, #Nominal Aitchison and Aitken kernel
    u_aitken = function(A, B, bw, C) {sum(ifelse(A == B, 1, bw))}, #Nominal Aitken kernel
    o_wangvanryzin = function(A, B, bw, C) {sum(ifelse(all(A == B), 1 - bw, (1/2) * (1 - bw) * (bw^abs(A - B))))}, #Ordinal Wang & van Ryzin kernel
    o_aitchisonaitken = function(A, B, bw, C) {sum(ifelse(A == B, 1, bw^(abs(A - B))))}, #Ordinal Aitchison and Aitken kernel
    o_aitken = function(A, B, bw, C) {sum(ifelse(A == B, bw, (1 - bw) / (2^abs(A - B))))}, #Ordinal Aitken kernel
    o_habbema = function(A, B, bw, C) {sum(bw^((abs(A - B))^2))} #Ordinal Habbema kernel
  )
 kernel.distance.sum <- function(cIND, uIND, oIND, cFUN, uFUN, oFUN, A, B, bw, df) {
  get_function <- function(FUN_name, type) {
     if (!FUN_name %in% names(kernel_functions)) {
       stop(paste("Invalid", type, "function name:", FUN_name))
     }
     return(kernel_functions[[FUN_name]])
   }
  
  compute_distance <- function(FUN, ind) {
    if (length(ind) == 0) return(0)
    sum(FUN(A[ind], A[ind], bw[ind], df[, ind]) + FUN(B[ind], B[ind], bw[ind], df[, ind]) - FUN(A[ind], B[ind], bw[ind], df[, ind]) - FUN(B[ind], A[ind], bw[ind], df[, ind]))
  }
  d_c <- if (cIND > 0) compute_distance(get_function(cFUN, "cFUN"), 1:cIND) else 0
  d_u <- if (uIND > 0) compute_distance(get_function(uFUN, "uFUN"), (cIND + 1):uIND) else 0
  d_o <- if (oIND > 0) compute_distance(get_function(oFUN, "oFUN"), (cIND + uIND + 1):oIND) else 0
  return(d_c + d_u + d_o)
}
  
  for (i in 1:ncol(combinations)) {
    row1 <- combinations[1, i]
    row2 <- combinations[2, i]
    distance <- kernel.distance.sum(cIND, uIND, oIND, cFUN, uFUN, oFUN, df[row1, ], df[row2, ], bw, df)
    distances[row1, row2] <- as.numeric(distance)
    distances[row2, row1] <- as.numeric(distance)
  }
  if (stan == FALSE) return(distances)
  if (stan == TRUE) {
    min <- min(distances)
    max <- max(distances)
    standardized <- (distances - min) / (max - min)
    return(standardized)
  }
}
