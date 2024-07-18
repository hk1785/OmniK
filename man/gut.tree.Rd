\name{gut.tree}
\alias{gut.tree}
\docType{data}
\title{
A phylogenetic tree for the gut microbiome data (Yanai et al., 2024, Nat Commun)
}
\description{
A phylogenetic tree for the gut microbiome data (Yanai et al., 2024, Nat Commun)
}
\usage{data("gut.tree")}
\format{
  The format is:
List of 5
 $ edge       : int [1:836, 1:2] 420 421 421 422 422 420 423 423 424 425 ...
 $ edge.length: num [1:836] 0.0383 0.0789 0.0802 0.0451 0.0477 ...
 $ Nnode      : int 418
 $ node.label : chr [1:418] "" "0.752" "0.992" "0.238" ...
 $ tip.label  : chr [1:419] "e920f96de1114d11639e145ceaddfaa3" "160c129cec3f44126b24eed8dd8c4621" "087a0803b5a078c7c547580b445f0903" "0c0603a2fdeae1308fceaaf717df0de2" ...
 - attr(*, "class")= chr "phylo"
 - attr(*, "order")= chr "cladewise"
}
\references{
Yanai H, Park B, Koh H. Jang HJ, Vaughan KL, Tanaka-Yano M, Aon M, Blanton M, Messaoudi I, Diaz-Ruiz A, Mattison JA, Beerman I. Short-term periodic restricted feeding elicits metabolome-microbiome signatures with sex dimorphic persistence in primate intervention. Nat Commun, 15(1088), 2024

Koh H. A general kernel machine regression framework using principal component analysis for jointly testing main and interaction effects. (In Review)
}
\examples{
data(gut.tree)
}
\keyword{datasets}
