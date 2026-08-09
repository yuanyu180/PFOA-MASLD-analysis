library(limma)
library(dplyr)
library(ggplot2)
data <- read.table("GSE66676.normalize.txt", header = TRUE, row.names = 1, sep = "\t")
dim(data)
head(data[, 1:5])
sample_names <- colnames(data)
group <- ifelse(grepl("_Control", sample_names), "Control", "Treat")
group <- factor(group, levels = c("Control", "Treat"))
table(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(data, design)
contrast.matrix <- makeContrasts(Treat_vs_Control = Treat - Control,
                                 levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, number = Inf, adjust.method = "BH")
head(results)
results$diffexpressed <- "NO"
results$diffexpressed[results$logFC > 1 & results$adj.P.Val < 0.05] <- "UP"
results$diffexpressed[results$logFC < -1 & results$adj.P.Val < 0.05] <- "DOWN"
table(results$diffexpressed)
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val),
                                    color = diffexpressed)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
print(volcano_plot)
write.csv(results, "differential_expression_results.csv", row.names = TRUE)
significant_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
write.csv(significant_genes, "significant_genes.csv", row.names = TRUE)
top_genes <- rownames(results[1:50, ])
heatmap_data <- as.matrix(data[top_genes, ])
heatmap_data <- t(scale(t(heatmap_data)))
library(pheatmap)
pheatmap(heatmap_data,
         annotation_col = data.frame(Group = group, row.names = colnames(data)),
         show_rownames = FALSE,
         main = "Top 50 Differentially Expressed Genes")
