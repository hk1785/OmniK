\name{D.to.K}
\alias{D.to.K}
\title{
Distance to kernel
}
\description{
This function converts pairwise (subject-by-subject, n by n) distance matrices to pairwise (subject-by-subject, n by n) kernel matrices.
}
\usage{
D.to.K(D)
}
\arguments{
  \item{D}{
A list of pairwise (subject-by-subject, n by n) distance matrices
}
}
\value{
A list of pairwise (subject-by-subject, n by n) kernel matrices
}
\references{
Koh H. A general kernel machine regression framework using principal component analysis for jointly testing main and interaction effects. (In Review)
}
\author{
Hyunwook Koh <hyunwook.koh@stonybrook.edu>
}

\examples{
##############################################################
# Example 1: Gut microbiome (Yanai et al., 2024, Nat Commun) #
##############################################################

gut.Ks <- lapply(gut.Ds, FUN = function(d) D.to.K(d))

####################################################################
# Example 2: Oral microbiome (Park et al., 2023, BMC Microbiology) #
####################################################################

oral.Ks <- lapply(oral.Ds, FUN = function(d) D.to.K(d))
}
