# Load your harmonized data
data <- read.table("results/harmonized.tsv", header=TRUE)

# Calculate lambda
lambda <- median(qchisq(data$P, df=1, lower.tail=FALSE)) / qchisq(0.5, df=1)

# Correct p-values
data$P_GC <- pchisq(qchisq(data$P, df=1, lower.tail=FALSE)/lambda, df=1, lower.tail=FALSE)

# Save corrected data
write.table(data, "results/harmonized_GC.tsv", sep="\t", row.names=FALSE)