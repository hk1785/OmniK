\name{OmniK}
\alias{OmniK}
\title{
OmniK
}
\description{
This function implements OmniK for the kernel machine regression analysis to jointly test main and interaction effects.
}
\usage{
OmniK(Y, Tr, Z = NULL, Ks, out.type = c("continuous", "binary"), n.perm = 3000)
}
\arguments{
  \item{Y}{
A numeric vector for responses.
}
  \item{Tr}{
A numeric vector for treatments.
}
  \item{Z}{
A matrix for covariates (e.g., age, sex).
}
  \item{Ks}{
A list of pairwise (subject-by-subject, n by n) kernel matrices.
}
  \item{out.type}{
The type of the response variable: "binary" or "continuous".
}
  \item{n.perm}{
The number of permutations (Default: 3000). 
}
}

\value{
pvals.input.kernels: P-values for each input kernel.
pval.endogenous.kernels: P-values for each endogenous kernel on the main effects (M), interaction effects (I) or both of them (B).
OmniK.pval: P-value for OmniK.
}

\references{
Koh H. A general kernel machine regression framework using principal component analysis for jointly testing main and interaction effects. (In Review)
}
\author{
Hyunwook Koh <hyunwook.koh@stonybrook.edu>
}
\note{
A large number of permutations is recommended for sufficient convergence. 
}

\examples{

##############################################################
# Example 1: Gut microbiome (Yanai et al., 2024, Nat Commun) #
##############################################################

gut.Ks <- lapply(gut.Ds, FUN = function(d) D.to.K(d))

Y <- gut.meta$body.weight
Tr <- gut.meta$diet
Z <- gut.meta[,c("age", "sex")] 

set.seed(521)

gut.out <- OmniK(Y = Y, Tr = Tr, Z = Z, Ks = gut.Ks, out.type = "continuous", n.perm = 3000) 

gut.out$pvals.input.kernels
gut.out$pval.endogenous.kernels
gut.out$OmniK.pval

####################################################################
# Example 2: Oral microbiome (Park et al., 2023, BMC Microbiology) #
####################################################################

oral.Ks <- lapply(oral.Ds, FUN = function(d) D.to.K(d))

Y <- oral.meta$gingival.inflammation
Tr <- oral.meta$ecig.use
Z <- oral.meta[,c("age", "sex")] 

set.seed(521)

oral.out <- OmniK(Y = Y, Tr = Tr, Z = Z, Ks = oral.Ks, out.type = "binary", n.perm = 3000) 

oral.out$pvals.input.kernels
oral.out$pval.endogenous.kernels
oral.out$OmniK.pval

}
