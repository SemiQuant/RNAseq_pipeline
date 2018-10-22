library("DESeq2")



countdata = apply(all_counts, 2, function(x) as.integer(x))
rownames(countdata)= rownames(all_counts)

#lets remove those with many zeros
remove = apply(countdata, 2, function(x) table(x == 0)[2] > 500)
remove = remove[remove == "TRUE"]

countdata = countdata[, !(colnames(countdata)%in%names(remove))]

#lets make condition same
mycoldata <- read.table("/Users/jdlim/Desktop/CASS_RNA_test/Status.txt", header = T)
mycoldata$condition = as.factor(mycoldata$condition)
mycoldata$file=gsub("_R1_001.fastq.gz.abundances.cxb", "", mycoldata$file)
mycoldata = mycoldata[mycoldata$file%in%colnames(countdata),]

countdata = countdata[,colnames(countdata)%in%mycoldata$file]

dim(countdata)
dim(mycoldata)

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
min(res$padj, na.rm = T)
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

