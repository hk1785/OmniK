# OmniK

Title: A general kernel machine regression framework using principal component analysis for jointly testing main and interaction effects

Version: 1.0

Date: 2024-07-18

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: This R package provides facilities for OmniK to conduct kernel machine regression analysis analysis for jointly testing main and interaction effects.

Depends: R(>= 4.1.1)

License: GPL-2

NeedsCompilation: no

URL: https://github.com/hk1785/OmniK

## Reference

* Koh H. A general kernel machine regression framework using principal component analysis for jointly testing main and interaction effects. (_In Review_)

## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/OmniK/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).

## Installation

```
library(devtools)
install_github("hk1785/OmniK", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------

## :mag: OmniK

### Description
This function implements OmniK for the kernel machine regression analysis to jointly test main and interaction effects.

### Usage
```
OmniK(Y, Tr, Z = NULL, Ks, out.type = c("continuous", "binary"), n.perm = 3000)
```

### Arguments
* _Y_ - A numeric vector for responses.
* _Tr_ - A numeric vector for treatments.
* _Z_ - A matrix for covariates (e.g., age, sex).
* _Ks_ - A list of pairwise (subject-by-subject, n by n) kernel matrices.
* _out.type_ - The type of the response variable: "binary" or "continuous".
* _n.perm_ - The number of permutations (Default: 3000). 

### Values
* _pvals.input.kernels_ - P-values for each input kernel.
* _pval.endogenous.kernels_ - P-values for each endogenous kernel on the main effects (M), interaction effects (I) or both of them (B).
* _OmniK.pval_ - P-value for OmniK.

### Example
Import OmniK
```
library('OmniK')
```
Example 1: Gut microbiome (Yanai et al., 2024, Nat Commun)
```
gut.Ks <- lapply(gut.Ds, FUN = function(d) D.to.K(d))

Y <- gut.meta$body.weight
Tr <- gut.meta$diet
Z <- gut.meta[,c("age", "sex")] 

set.seed(521)

gut.out <- OmniK(Y = Y, Tr = Tr, Z = Z, Ks = gut.Ks, out.type = "continuous", n.perm = 3000) 

gut.out$pvals.input.kernels
gut.out$pvals.endogenous.kernels
gut.out$OmniK.pval
```
Example 2: Oral microbiome (Park et al., 2023, BMC Microbiology)
```
oral.Ks <- lapply(oral.Ds, FUN = function(d) D.to.K(d))

Y <- oral.meta$gingival.inflammation
Tr <- oral.meta$ecig.use
Z <- oral.meta[,c("age", "sex")] 

set.seed(521)

oral.out <- OmniK(Y = Y, Tr = Tr, Z = Z, Ks = oral.Ks, out.type = "binary", n.perm = 3000) 

oral.out$pvals.input.kernels
oral.out$pvals.endogenous.kernels
oral.out$OmniK.pval
```
