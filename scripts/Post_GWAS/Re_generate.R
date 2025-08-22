#!/usr/bin/env Rscript

#install packages
#install.packages("qqman")
#install.packages("data.table")

#load packages
library(qqman)
library(data.table)

# Load harmonized data
gwas <- fread("results/harmonized_GC.tsv")

# 1. Manhattan Plot (All Chromosomes)
png("results/re-plots/manhattan.png", width=1200, height=400)
manhattan(gwas, chr="CHR", bp="POS", p="P_GC", snp="SNP",
          main="Breast Cancer GWAS - Harmonized Variants",
          suggestiveline=-log10(1e-5),
          genomewideline=-log10(5e-8))
dev.off()

# 2. QQ Plot (Check Inflation)
lambda <- median(qchisq(gwas$P_GC, df=1, lower.tail=F)) / qchisq(0.5, df=1)
png("results/re-plots/qqplot.png", width=600, height=600)
qq(gwas$P_GC, main=paste("QQ Plot (λ =", round(lambda, 3), ")"),
   col="blue", cex=1.5)
grid()
dev.off()

# 3. Effect Size Distribution
png("results/re-plots/effect_sizes.png", width=800, height=600)
hist(gwas$BETA, breaks=30, col="#4E79A7",
     main="Distribution of Effect Sizes",
     xlab="Beta (Effect Size)")
dev.off()