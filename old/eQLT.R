# eQLT analysis
# http://www.bioconductor.org/packages/release/bioc/vignettes/iBMQ/inst/doc/iBMQ.pdf
require(iBMQ)


data(snp)
data(gene)

# This function computes the MCMC algorithm to produce Posterior Probabilities of Association (PPA) for eQTL mapping.
# preferably in the order of 100,000, with a burn-in of 50,000)
PPA <- eqtlMcmc(snp, gene, n.iter=100,burn.in=100,n.sweep=20,mc.cores=6, RIS=TRUE)

snp_data <- snp@assayData$call
snp_data <- new('SnpSet', call = snp_data)
# call (a matrix of genotypic calls, with features (SNPs) corresponding to rows and samples to columns)
#  experimentData = [MIAME], phenoData = [AnnotatedDataFrame], annotation = [character], protocolData = [AnnotatedDataFrame],  call = [matrix], callProbability = [matrix], ...)
gene_data <- gene@assayData$exprs
gene_data <- ExpressionSet(assayData=gene_data)
#trying with just data frame
PPA <- eqtlMcmc(snp_data, gene_data, n.iter=100, burn.in=100, n.sweep=20, mc.cores=6, RIS=TRUE)



# a FDR of 0.1 (corresponding to 10%)
cutoff <- calculateThreshold(PPA, 0.01)

# the PPA optimal cutoff == 0.74. This means that all the eQTLs with a PPA above 0.74 are significant eQTLs
eqtl <- eqtlFinder(PPA, cutoff)

# Classifying the eQTLs
data(snppos)
data(genepos)
# We need to specify a cutoff value(in base pair) corresponding to the threshold where a eQTL
# is considered a cis-eQTL. In this case, we will use a cutoff of 1 MB (1000000 pb).
eqtltype <- eqtlClassifier(eqtl, snppos, genepos, 1000000)

# Visualizing the result
library(ggplot2)
ggplot(eqtltype, aes(y=GeneStart, x=MarkerPosition)) +
  geom_point(aes(y=GeneStart, x=MarkerPosition, color = PPA), size = 1.5)+
  facet_grid(GeneChrm~MarkerChrm)+theme_bw(base_size = 12, base_family = "")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())+scale_x_reverse()

# Finding eQTL hotspots
hotspot <- hotspotFinder(eqtltype, 10)











