library(bio3d)

mypdbfile <- "path_to_pdb"
mydcdfile <- "path_to_dcd"

dcd <- read.dcd(mydcdfile)
pdb <- read.pdb(mypdbfile)

ca.inds <- atom.select(pdb, elety="CA")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)

pc <- pca.xyz(xyz[,ca.inds$xyz])

hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)
