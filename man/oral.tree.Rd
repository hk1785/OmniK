\name{oral.tree}
\alias{oral.tree}
\docType{data}
\title{
A phylogenetic tree for the oral microbiome data (Park et al., 2023, BMC Microbiology)
}
\description{
A phylogenetic tree for the oral microbiome data (Park et al., 2023, BMC Microbiology)
}
\usage{data("oral.tree")}
\format{
  The format is:
List of 5
 $ edge       : int [1:3013, 1:2] 1548 1549 1550 1550 1549 1551 1552 1553 1554 1555 ...
 $ edge.length: num [1:3013] 3.49e-02 1.81e-01 5.00e-09 2.37e-03 3.27e-02 ...
 $ Nnode      : int 1467
 $ node.label : chr [1:1467] "0.588" "0.985" "0.591" "0.961" ...
 $ tip.label  : chr [1:1547] "376f2df812b1f4085a8ebceea456e562" "03b349368bcce56cde78611b3b2c4eb7" "bc0437f93952da3a2ac9c9624c471ce9" "d41e476312d2c73de292b0885fd2a3bb" ...
 - attr(*, "class")= chr "phylo"
 - attr(*, "order")= chr "cladewise"
}
\references{
Park B, Koh H, Patatanian M, Reyes-Caballero H, Zhao N, Meinert J, Holbrook JT, Leinbach LI, Biswal S. The mediating roles of the oral microbiome in saliva and subgingival sites between e-cigarette smoking and gingival inflammation. BMC Microbiol, 23(35), 2023

Koh H. A general kernel machine regression framework using principal component analysis for jointly testing main and interaction effects. (In Review)
}
\examples{
data(oral.tree)
}
\keyword{datasets}
