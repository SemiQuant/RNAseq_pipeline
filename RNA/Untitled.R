#RNA analysis

all <- read.table("/Users/jdlim/Library/Mobile Documents/com~apple~CloudDocs/Work/RNA/RNA_optimization experiment/RNA_sheet1.txt",
                  sep = "\t", header = T)

colnames(all)
#
# require(dplyr)
#
# cdata <- ddply(all, c("Arm", "Bio.replicate"), summarise,
#                mean = mean(change, na.rm=TRUE))
#
#


all.means <- aggregate(all,
                         by = all[c("Arm", "Type", "Bio.replicate")], FUN=median)

colnames(all.means)


all.medians <- aggregate(all[c("X18s", "GAPDH", "X16s", "SigA_in.run")],
                         by = all[c("Arm", "Type", "Bio.replicate")], FUN=median)



all.var <- aggregate(all,
                       by = all[c("Arm", "Type", "Bio.replicate")], FUN=var)

a = all.var[c("Arm", "Type", "Bio.replicate",
              "X18s..median.conc.original.", "X18s_cDNA_neg..median.conc.original.", "GAPDH..median.conc.original.",
              "GAPDH_cDNA_neg..median.conc.original.", "X16s..median.conc.original.", "X16s_cDNA_neg..median.conc.original.",
              "SigA..median.conc.original.", "SigA_cDNA_neg..median.conc.original." )]


b= all.means[c("Arm", "Type", "Bio.replicate",
             "X18s..median.conc.original.", "X18s_cDNA_neg..median.conc.original.", "GAPDH..median.conc.original.",
             "GAPDH_cDNA_neg..median.conc.original.", "X16s..median.conc.original.", "X16s_cDNA_neg..median.conc.original.",
             "SigA..median.conc.original.", "SigA_cDNA_neg..median.conc.original." )]



write.table(all.means, "/Users/jdlim/Desktop/temp_means_fixedVal", sep = "\t", row.names = F)




#
#
#
#
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("EasyqpcR")
# require(EasyqpcR)
#
#
# # biocLite("HTqPCR")
# # require(HTqPCR)
# # path = "/Users/jdlim/Library/Mobile Documents/com~apple~CloudDocs/Work/RNA/RNA_optimization experiment/lightCyclerCalcs/backup/rawruns/rRNAs_cycles.txt"
# # test <- read.delim(path)
# # raw <- readCtData(path, n.features = 384, format = "LightCycler")
#
#
# require(ReadqPCR)
# cycData <- read.LC480(file = path)
#
#

