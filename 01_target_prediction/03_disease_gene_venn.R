library(ggvenn)
files <- list.files(pattern = "*.txt$")
geneList <- list()
for (inputFile in files) {
  rt <- read.table(inputFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  fileName <- gsub("\\.txt$", "", inputFile)
  geneList[[fileName]] <- unique(rt[, 1])
}
union_genes <- Reduce(union, geneList)
intersection_genes <- Reduce(intersect, geneList)
write.table(union_genes, file = "union_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(intersection_genes, file = "intersection_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
color_count <- length(geneList)
fill_colors <- c("#E41A1C", "#1E90FF", "#FF8C00", "#32CD32", "#9400D3", "#A52A2A")[1:color_count]
set_name_color <- rep("black", color_count)
pdf("venn.pdf", width = 6, height = 6)
ggvenn(
  geneList,
  show_percentage = TRUE,
  stroke_color = "white",
  stroke_size = 0.5,
  fill_color = fill_colors,
  set_name_color = set_name_color,
  set_name_size = 7,
  text_size = 4.5
)
dev.off()
