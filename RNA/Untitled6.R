file <- "/Users/jdlim/Desktop/IS6110/CASS_isolates/IS_results/k/k_compiled_table_IS6110_edited.txt"

table <- read.table(file, header=T)
rownames(table) <- table$isolate
table <- table[,-1]

count_question <- apply(table, 1, function(x) sum(x == "?"))
count_all <- apply(table, 1, function(x) sum(x == "?" | x == "+" | x == "*"))


table$count_all <- count_all
table$count_question <- count_question
write.table(table, file)





#check fit
a <- read.table("/Users/jdlim/Desktop/temp.txt", header = T)

unique(a$Lineage)
a <- a[a$Lineage == "lineage2.2.1", ]
plot(a$Number_of_bands_RFLP, (a$NumberofIS6110elementsWGS_CDC_genome_standar_settings-a$NumberofIS6110elementsWGS_CDC_genome_standar_settings_lowQual))
abline(lm((a$NumberofIS6110elementsWGS_CDC_genome_standar_settings-a$NumberofIS6110elementsWGS_CDC_genome_standar_settings_lowQual) ~a$Number_of_bands_RFLP  ))
p <- lm((a$NumberofIS6110elementsWGS_CDC_genome_standar_settings-a$NumberofIS6110elementsWGS_CDC_genome_standar_settings_lowQual)~a$Number_of_bands_RFLP)
anova(p)[1,'Pr(>F)']
anova(p)



#make gel
#will need to edit values as gel does not seperate at even speed
#can calculate from the lab ones

table <- read.table("/Users/jdlim/Desktop/IS6110/CASS_isolates/IS_results/H37/iteration1/test/compiled_table_IS6110_counts.txt", header=T)
rownames(table) <- table$isolate
table <- table[,-1]
table <- table[,-ncol(table)]

positions <- colnames(table)
positions <- gsub("X", "", positions)
positions <- sub("\\..*","", positions)

table_test <- t(table)
row.names(table_test) <- positions
table_test[table_test == "+"] <- 1
table_test[table_test == "*"] <- 1
table_test[table_test == "?"] <- 0
table_test[table_test == "-"] <- 0

max <- 4411709+1000
max <- log2(sort(max))#*log10(sort(max))

plot(x = -10, y = -10, xlim=c(0,ncol(table_test)), ylim=c(0, max), pch = NULL, xlab = "Isolate", ylab = "", xaxt='n' )
axis(1, labels = colnames(table_test), at = seq(1, ncol(table_test)), cex.axis=0.2, las = 2)

for (i in 1:ncol(table_test)){

  temp <- table_test[,c(i)]
  temp <- temp[temp==1]
  temp <- as.numeric(names(temp))

  temp3 <- NULL

  for ( position in 1:length(temp)-1){
    temp2 <- temp[position+1] - temp[position]
    temp3 <- append(temp3, temp2)
  }

  temp4 <- 4411709 - temp[length(temp)] + temp[1]
  temp3 <- append(temp3, temp4)

  temp3 <- log2(sort(temp3))#*log10(sort(temp3))
  temp3 <- max - temp3

  par(new=T)
  plot(y = temp3, x = rep(i, length(temp3)), ylim = c(0, max), xlim=c(0,ncol(table_test)), pch = "-", xlab = "", ylab = "", axes = F)
}






