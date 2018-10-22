a=read.table("/Users/jdlim/Desktop/CASS_RNA_test/testing/gene_exp.diff", header = T)
b=a[a$gene=="ansA",]





c = read.table("/Users/jdlim/Desktop/CASS_RNA_test/testing/genes.fpkm_tracking", header = T)
d = c[, grepl(pattern = "FPKM", x = colnames(c)) | colnames(c) == "gene_short_name"]
row_sub = apply(d[c(2:ncol(d))], 1, function(row) all(row !=0 ))
##Subset as usual
d = d[row_sub,]



require(cummeRbund)
cuff<-readCufflinks(dir = "/Users/jdlim/Desktop/CASS_RNA_test/testing/")

disp<-dispersionPlot(genes(cuff))
disp
# To assess the distributions of FPKM scores across samples, you can use the csDensity plot (Figure 1).
dens<-csDensity(genes(cuff))
dens
densRep<-csDensity(genes(cuff),replicates=T)
densRep
#
#
# s<-csScatterMatrix(genes(cuff))
# s




a <- read.table("/Users/jdlim/Desktop/CASS_RNA_test/cuffnorm_test/genes.count_table", header = T)
row_sub <- apply(a[c(2:ncol(a))], 1, function(row) all(row !=0 ))
##Subset as usual
b <- a[row_sub,]

c=colSums(a[,-1])


boxplot(a[,-1])


#######################################################################
#### Select top 200 features using Median Absolute Deviation (MAD) ####
#######################################################################
mads = apply(a[,-1],1,mad)
# Plot histogram of MAD before feature selection
# pdf('histogram_median_absolute_deviations_all_genes.pdf')
hist(mads,xlab='Median Absolute Deviation (MAD)')
# dev.off()
d1 = a[,-1][order(mads,decreasing=T)[1:200],]

# Plot histograms of MAD after feature selection
# pdf('histogram_median_absolute_deviations_selected_genes.pdf')
hist(mads[order(mads,decreasing=T)[1:200]],xlab='Median Absolute Deviation (MAD)')
# dev.off()




















countfile <- "/Users/jdlim/Desktop/CASS_RNA_test/cuffnorm_test/genes.fpkm_table"
countdata <- read.table(countfile, header = TRUE, row.names = 1, stringsAsFactors = F)
library("DESeq2")

mycoldata <- read.table("/Users/jdlim/Desktop/CASS_RNA_test/Status.txt", header = T)
mycoldata$condition = as.factor(mycoldata$condition)

countdata = apply(countdata, 2, function(x) as.integer(x))

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = mycoldata, design = ~condition)


dds <- DESeq(dds)
# Now, let???s use the results() function to pull out the results from the dds object. Let???s
# re-order by the adjusted p-value.
# Get differential expression results
res <- results(dds)
# head(res)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# head(res)
# Combine DEseq results with the original counts data. Write significant results to a file.
sig <- subset(res, padj < 0.05)
sig

# Transform
rld <- rlogTransformation(dds)
# Principal components analysis
plotPCA(rld, intgroup = "condition")
# Hierarchical clustering analysis let's get the actual values for the first
# few genes
# head(assay(rld))
## now transpose those
# t(head(assay(rld)))
## now get the sample distances from the transpose of the whole thing
# dist(t(assay(rld)))
sampledist <- dist(t(assay(rld)))
plot(hclust(sampledist))
# Let???s plot a heatmap.
# ?heatmap for help
# sampledist
# as.matrix(sampledist)
sampledistmat <- as.matrix(sampledist)

temp = data.frame(colnames(sampledistmat))
temp$order = seq(1, nrow(temp))
temp = merge(temp, mycoldata[c(1,3)], by.x="colnames.sampledistmat.", by.y="sample_id")
temp=temp[order(temp$order),]
temp$jointName = paste(temp$condition, temp$colnames.sampledistmat.)
colnames(sampledistmat) = temp$jointName; rownames(sampledistmat) = temp$jointName


heatmap(sampledistmat, col = heat.colors(256))
# col = colorRampPalette(c( "blue", "red"))(n = 100) )


library(gplots)
heatmap.2(sampledistmat, col=heat.colors(30),
          key=T, keysize=1.5,
          trace="none")






#heatmap of counts
countfile <- "/Users/jdlim/Desktop/CASS_RNA_test/cuffnorm_test/genes.count_table"
countdata <- read.table(countfile, header = TRUE, row.names = 1, stringsAsFactors = F)
summary(countdata)


countdata_sigA <- countdata[rownames(countdata) == "Rv2703",]
summary(countdata_sigA)

library(gplots)
heatmap.2(as.matrix(countdata), col=heat.colors(100),
          key=T, keysize=1.5,
          trace="none")






require(cummeRbund)
cuff<-readCufflinks(dir = "/Users/jdlim/Desktop/CASS_RNA_test/testing/")
myGeneIds <- unlist(read.table("/Users/jdlim/bioinfomatics/ecamber/datasets/mtu2/anns_parsed/refseq/H37Rv.txt", sep = "\t")[1])

myGenes<-getGenes(cuff,myGeneIds)

csHeatmap(isoforms(myGenes),cluster="both",labRow=F, pseudocount = T)


