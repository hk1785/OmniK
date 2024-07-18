OmniK <-
function(Y, Tr, Z = NULL, Ks, out.type = c("continuous", "binary"), n.perm = 3000) {
  
  H <- length(Ks)

  if(!is.null(Z) ) {
    meta.d <- as.data.frame(cbind(Y, Tr, Z))
    ind.na <- apply(meta.d, 1, function(x) sum(is.na(x)))
    meta.d <- as.data.frame(model.matrix(~ ., data = meta.d)[,-1])
    Y <- as.numeric(meta.d[, 1])
    Tr <- as.numeric(meta.d[, 2])
    Z <- as.matrix(meta.d[,-c(1, 2)])
    for (i in 1:H) {
      Ks[[i]] <- Ks[[i]][ind.na == 0, ind.na == 0]
    }
  }
  
  if(is.null(Z) ) {
    meta.d <- as.data.frame(cbind(Y, Tr))
    ind.na <- apply(meta.d, 1, function(x) sum(is.na(x)))
    meta.d <- as.data.frame(model.matrix(~ ., data = meta.d)[,-1])
    Y <- as.numeric(meta.d[, 1])
    Tr <- as.numeric(meta.d[, 2])
    for (i in 1:H) {
      Ks[[i]] <- Ks[[i]][ind.na == 0, ind.na == 0]
    }
  }
  
  n <- nrow(Ks[[1]])
  
  ran.perm.ind <- list()
  for (i in 1:n.perm) {
    ran.perm.ind[[i]] <- sample(1:n, replace = FALSE)
  }
  
  Ms <- lapply(Ks, function(x) {
    svd.out <- svd(x)
    ind <- which(svd.out$d > 0)
    svd.out$u[,ind] %*% sqrt(diag(svd.out$d[ind]))})
  
  if(is.null(Z) ) {
    if (out.type == "continuous") {
      Y <- Y - glm(Y ~ Tr, family = gaussian())$fitted.values
    }
    
    if (out.type == "binary") {
      Y <- Y - glm(Y ~ Tr, family = binomial())$fitted.values
    }
  }
  
  if(!is.null(Z) ) {
    if (out.type == "continuous") {
      Y <- Y - glm(Y ~ Tr + Z, family = gaussian())$fitted.values
    }
    
    if (out.type == "binary") {
      Y <- Y - glm(Y ~ Tr + Z, family = binomial())$fitted.values
    }
  }
  
  # Main
  
  Ks.M <- lapply(Ms, function(x) x %*% t(x))
  
  T.M.obs.list <- sapply(Ks.M, function(x) as.numeric(t(Y) %*% x %*% Y))
  
  T.M.perm.list <- list()
  for (i in 1:H) {
    T.M.perm.list[[i]] <- sapply(ran.perm.ind, function(x) as.numeric(t(Y[x]) %*% Ks.M[[i]] %*% Y[x])) 
  }    
  
  itembyitem.M.pvals <- numeric()
  for (i in 1:H) {
    itembyitem.M.pvals[i] <- sum(T.M.perm.list[[i]] > T.M.obs.list[i])/(n.perm + 1)
  }
  
  # Interaction
  
  Ks.I <- lapply(Ms, function(x) (x*Tr) %*% t(x*Tr))
  
  T.I.obs.list <- sapply(Ks.I, function(x) as.numeric(t(Y) %*% x %*% Y))
  
  T.I.perm.list <- list()
  for (i in 1:H) {
    T.I.perm.list[[i]] <- sapply(ran.perm.ind, function(x) as.numeric(t(Y[x]) %*% Ks.I[[i]] %*% Y[x])) 
  }    
  
  itembyitem.I.pvals <- numeric()
  for (i in 1:H) {
    itembyitem.I.pvals[i] <- sum(T.I.perm.list[[i]] > T.I.obs.list[i])/(n.perm + 1)
  }
  
  # Main and Interaction
  
  Ks.B <- lapply(Ms, function(x) cbind(x, (x*Tr)) %*% t(cbind(x, (x*Tr))))
  
  T.B.obs.list <- sapply(Ks.B, function(x) as.numeric(t(Y) %*% x %*% Y))
  
  T.B.perm.list <- list()
  for (i in 1:H) {
    T.B.perm.list[[i]] <- sapply(ran.perm.ind, function(x) as.numeric(t(Y[x]) %*% Ks.B[[i]] %*% Y[x])) 
  }    
  
  itembyitem.B.pvals <- numeric()
  for (i in 1:H) {
    itembyitem.B.pvals[i] <- sum(T.B.perm.list[[i]] > T.B.obs.list[i])/(n.perm + 1)
  }
  
  # P-values for input kernels
  
  itembyitem.pvals <- cbind(itembyitem.M.pvals, itembyitem.I.pvals, itembyitem.B.pvals)
  MinP.itembyitem.obs <- apply(itembyitem.pvals, 1, min)
  
  MinP.itembyitem.perm.list <- list()
  for (j in 1:H) {
    T.X.perm.list <- list(T.M.perm.list[[j]], T.I.perm.list[[j]], T.B.perm.list[[j]])
    MinP.itembyitem.perm <- numeric()
    for (i in 1:n.perm) {
      MinP.itembyitem.perm[i] <- min(sapply(T.X.perm.list, function(v) sum(abs(v[-i]) > abs(v[i]))/n.perm))
    }  
    MinP.itembyitem.perm.list[[j]] <- MinP.itembyitem.perm
  }
  
  MinP.itembyitem.pval <- numeric()
  for (j in 1:H) {
    MinP.itembyitem.pval[j] <- (sum(MinP.itembyitem.perm.list[[j]] < MinP.itembyitem.obs[j]) + 1)/(n.perm + 1)
  }
  
  # P-values for endogenous kernels
  
  MinP.OmniK.M.obs <- min(itembyitem.M.pvals)
  
  MinP.OmniK.M.perm <- numeric()
  for (i in 1:n.perm) {
    MinP.OmniK.M.perm[i] <- min(sapply(T.M.perm.list, function(v) sum(abs(v[-i]) > abs(v[i]))/n.perm))
  }  
  
  MinP.OmniK.M.pval <- (sum(MinP.OmniK.M.perm < MinP.OmniK.M.obs) + 1)/(n.perm + 1)
  
  MinP.OmniK.I.obs <- min(itembyitem.I.pvals)
  
  MinP.OmniK.I.perm <- numeric()
  for (i in 1:n.perm) {
    MinP.OmniK.I.perm[i] <- min(sapply(T.I.perm.list, function(v) sum(abs(v[-i]) > abs(v[i]))/n.perm))
  }  
  
  MinP.OmniK.I.pval <- (sum(MinP.OmniK.I.perm < MinP.OmniK.I.obs) + 1)/(n.perm + 1)
  
  MinP.OmniK.B.obs <- min(itembyitem.B.pvals)
  
  MinP.OmniK.B.perm <- numeric()
  for (i in 1:n.perm) {
    MinP.OmniK.B.perm[i] <- min(sapply(T.B.perm.list, function(v) sum(abs(v[-i]) > abs(v[i]))/n.perm))
  }  
  
  MinP.OmniK.B.pval <- (sum(MinP.OmniK.B.perm < MinP.OmniK.B.obs) + 1)/(n.perm + 1)
  
  # OmniK
  
  MinP.OmniK.obs <- min(itembyitem.pvals)
  
  T.X.perm.list <- c(T.M.perm.list, T.I.perm.list, T.B.perm.list)
  MinP.OmniK.perm <- numeric()
  for (i in 1:n.perm) {
    MinP.OmniK.perm[i] <- min(sapply(T.X.perm.list, function(v) sum(abs(v[-i]) > abs(v[i]))/n.perm))
  }
  
  MinP.OmniK.pval <-  (sum(MinP.OmniK.perm < MinP.OmniK.obs) + 1)/(n.perm + 1)
  
  pvals.endogenous.kernels <- c(MinP.OmniK.M.pval, MinP.OmniK.I.pval, MinP.OmniK.B.pval)
  names(pvals.endogenous.kernels) <- c("M", "I", "B")
  
  output <- list(pvals.input.kernels = MinP.itembyitem.pval, 
                 pvals.endogenous.kernels = pvals.endogenous.kernels, 
                 OmniK.pval = MinP.OmniK.pval)
  
  return(output)
  
}
