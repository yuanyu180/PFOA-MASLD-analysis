Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
setwd(".")
library(Seurat)
library(dplyr)
library(scDblFinder)
library(future)
library(furrr)
library(harmony)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
samples <- c(
  "Control.GSM4041156",
  "Control.GSM4041159",
  "Control.GSM4041160",
  "Control.GSM4041151",
  "Control.GSM4041154",
  "Control.GSM5325534/",
  "Control.GSM5325535/",
  "Treat.GSM4041162",
  "Treat.GSM4041165",
  "Treat.GSM4041167",
  "Treat.GSM4041168",
  "Treat.GSM4041169",
  "Treat.GSM5325536/",
  "Treat.GSM5325537/"
)
sample_names <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6","Sample7", "Sample8", "Sample9", "Sample10", "Sample11", "Sample12", "Sample13", "Sample14")
process_sample <- function(sample_dir, sample_name) {
  message("Processing sample: ", sample_name)
  data <- Read10X(sample_dir)
  seu <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  seu$sample_id <- sample_name
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu[["log10GenesPerUMI"]] <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
  seu <- subset(seu, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &
                  percent.mt < 20 & nCount_RNA > 1000 & log10GenesPerUMI > 0.8)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce)
  seu$doublet.class <- ifelse(sce$scDblFinder.class == "singlet", "Singlet", "Doublet")
  seu$doublet.score <- sce$scDblFinder.score
  seu_filtered <- subset(seu, subset = doublet.class == "Singlet" & doublet.score < 0.8)
  seu_filtered <- RenameCells(seu_filtered, add.cell.id = sample_name)
  seu_filtered <- DietSeurat(seu_filtered, counts = TRUE, data = TRUE, assays = "RNA")
  return(seu_filtered)
}
message("Processing all samples...")
filtered_samples <- future_map2(samples, sample_names, process_sample, .options = furrr_options(seed = TRUE))
message("Merging samples...")
conbined <- Reduce(function(x, y) merge(x, y), filtered_samples)
message("Adding group information...")
conbined$group <- ifelse(conbined$sample_id %in% c("Sample1", "Sample2", "Sample3","Sample4", "Sample5", "Sample6", "Sample7"), "Healthy", "NAFLD")
conbined$condition <- conbined$group
cat("Group distribution:\n")
print(table(conbined$sample_id, conbined$group))
message("Processing cell cycle effects...")
cc.genes <- cc.genes.updated.2019
s.genes <- intersect(cc.genes$s.genes, rownames(conbined))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(conbined))
conbined <- CellCycleScoring(conbined,
                             s.features = s.genes,
                             g2m.features = g2m.genes,
                             set.ident = FALSE)
message("Data normalization and integration...")
message("Visualizing quality control metrics...")
p1 <- ggplot(conbined@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = group)) +
  geom_density_2d() +
  scale_color_manual(values = c("Healthy" = "blue", "NAFLD" = "red")) +
  labs(title = "Cell Count vs Feature Count", x = "Number of RNA counts", y = "Number of features (genes)") +
  theme_minimal()
p2 <- ggplot(conbined@meta.data, aes(x = percent.mt, fill = group)) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Healthy" = "blue", "NAFLD" = "red")) +
  labs(title = "Mitochondrial Gene Percentage Distribution", x = "Percentage of mitochondrial genes", y = "Cell Count") +
  theme_minimal()
p3 <- ggplot(conbined@meta.data, aes(x = log10GenesPerUMI, fill = group)) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Healthy" = "blue", "NAFLD" = "red")) +
  labs(title = "Log10 Genes per UMI Distribution", x = "Log10 Genes per UMI", y = "Cell Count") +
  theme_minimal()
ggsave("cell_count_vs_feature_count.png", plot = p1, width = 8, height = 6)
ggsave("mitochondrial_percentage_distribution.png", plot = p2, width = 8, height = 6)
ggsave("log10_genes_per_umi_distribution.png", plot = p3, width = 8, height = 6)
conbined_plot <- p1 + p2 + p3 + plot_layout(ncol = 1)
ggsave("quality_control_plots.png", plot = conbined_plot, width = 8, height = 12)
conbined <- NormalizeData(conbined)
conbined <- FindVariableFeatures(conbined, selection.method = "vst", nfeatures = 2000)
conbined <- ScaleData(conbined,
                      vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                      features = rownames(conbined))
conbined <- RunPCA(conbined, npcs = 30)
message("Running Harmony integration...")
conbined <- RunHarmony(conbined, group.by.vars = "sample_id", dims.use = 1:30)
conbined <- RunUMAP(conbined, reduction = "harmony", dims = 1:20)
conbined <- FindNeighbors(conbined, reduction = "harmony", dims = 1:20)
conbined <- FindClusters(conbined, resolution = 0.8)
message("Annotating cell types...")
umap_plot <- DimPlot(conbined,
                     reduction = "umap",
                     group.by = "seurat_clusters",
                     label = TRUE,
                     label.size = 6,
                     repel = TRUE,
                     pt.size = 1) +
  ggtitle("UMAP visualization of cell clusters") +
  theme(plot.title = element_text(hjust = 0.5))
print(umap_plot)
sample_umap <- DimPlot(conbined,
                       reduction = "umap",
                       group.by = "sample_id",
                       label = FALSE,
                       pt.size = 1) +
  ggtitle("UMAP visualization by sample") +
  theme(plot.title = element_text(hjust = 0.5))
umap_plot + sample_umap
table(Idents(conbined))
marker_genes <- list(
  "Monocytes" = c("CD14", "LYZ", "FCN1", "S100A9"),
  "B Cells" = c("CD79A", "MS4A1", "CD19", "CD79B","LRRC26"),
  "Endothelial Cells" = c("PECAM1", "VWF", "CD34", "CLDN5"),
  "NK Cells" = c("NKG7", "GNLY", "KLRD1", "NCAM1"),
  "T Cells" = c("CD3D", "CD3E", "CD2", "TRAC"),
  "Hepatocytes" = c("ALB", "APOA1", "HAMP", "TF"),
  "Stellate Cells" = c("ACTA2", "COL1A1", "PDGFRB", "DES"),
  "Neutrophils" = c("CSF3R", "FCGR3B", "S100A8", "S100A9"),
  "Macrophages" = c("CD68", "CD163", "MRC1", "C1QA","LILRA4"),
  "Kupffer Cells" = c("CD68", "MARCO", "VSIG4", "TIMD4"),
  "Dendritic Cells" = c("CD1C", "CLEC10A", "FCER1A", "CD83"),
  "Fibroblasts" = c("COL1A1", "DCN", "LUM", "PDGFRA"),
  "Cholangiocytes" = c("KRT19", "KRT7", "EPCAM", "SOX9")
)
p_dot <- DotPlot(conbined, features = unique(unlist(marker_genes)),
                 cluster.idents = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Gene Expression for Cell Type Annotation")
print(p_dot)
conbined=JoinLayers(conbined)
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("MAST")
library(MAST)
table(conbined@meta.data$seurat_clusters)
sc.s=sample(colnames(conbined),)
sc.s=subset(conbined,downsample=100)
table(sc.s@meta.data$seurat_clusters)
markers <- FindAllMarkers(sc.s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
current_clusters <- levels(conbined)
print("Current clusters:")
print(current_clusters)
new.cluster.ids <- c(
  "0" = "T Cells",
  "1" = "NK Cells",
  "2" = "Hepatocytes",
  "3" = "T Cells",
  "4" = "Endothelial Cells",
  "5" = "Endothelial Cells",
  "6" = "Monocytes",
  "7" = "Macrophages",
  "8" = "Monocytes",
  "9" = "Hepatocytes",
  "10" = "Hepatocytes",
  "11" = "Hepatocytes",
  "12" = "Cholangiocytes",
  "13" = "Endothelial Cells",
  "14" = "Hepatocytes",
  "15" = "B Cells",
  "16" = "Stellate Cells",
  "17" = "Stellate Cells" ,
  "18" = "Hepatocytes",
  "19" = "Hepatocytes" ,
  "20" = "Hepatocytes",
  "21" = "Hepatocytes",
  "22" = "NK Cells",
  "23" = "Hepatocytes",
  "24" = "Macrophages"
)
Idents(conbined) <- conbined$seurat_clusters
conbined <- RenameIdents(conbined, new.cluster.ids)
table(Idents(conbined))
conbined$celltype <- Idents(conbined)
umap_plot <- DimPlot(conbined,
                     reduction = "umap",
                     group.by = "celltype",
                     label = TRUE,
                     label.size = 6,
                     repel = TRUE,
                     pt.size = 1) +
  ggtitle("UMAP visualization of cell clusters") +
  theme(plot.title = element_text(hjust = 0.5))
print(umap_plot)
?saveRDS
saveRDS(conbined,file ="conbined.rds")
library(ggplot2)
marker_genes <- list(
  "Macrophages" = c("CD68", "CD163", "MRC1", "LILRA4"),
  "Monocytes" = c("CD14", "LYZ", "FCN1", "S100A9"),
  "Stellate Cells" = c("ACTA2", "COL1A1", "PDGFRB"),
  "B Cells" = c("CD79A", "MS4A1", "CD19", "CD79B"),
  "NK Cells" = c("NKG7", "GNLY", "KLRD1"),
  "T Cells" = c("CD3D", "CD3E", "CD2", "TRAC"),
  "Hepatocytes" = c("ALB", "APOA1", "TF"),
  "Cholangiocytes" = c("KRT19", "KRT7", "EPCAM", "SOX9"),
  "Endothelial Cells" = c("PECAM1", "VWF", "CD34", "CLDN5")
)
desired_order <- c('7', '24', '6', '8', '16', '15',  '1', '22', '0', '3', '11', '21', '14', '23', '19', '10', '2', '9', '17', '18','20', '12', '13', '4', '5')
conbined$seurat_clusters <- factor(conbined$seurat_clusters, levels = desired_order)
p_dot <- DotPlot(conbined, features = unique(unlist(marker_genes)),
                 group.by = "seurat_clusters") +
  ggtitle("Marker Gene Expression for Cell Type Annotation") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size(range = c(1, 6)) +
  ylab("Seurat Clusters") +
  xlab("Marker Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right")
print(p_dot)
p_dot <- DotPlot(conbined, features = unique(unlist(marker_genes)),
                 cluster.idents = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Gene Expression for Cell Type Annotation") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_dot)
selected_marker_genes <- c("BCL6","IL10","SHBG","NR4A2","CASP1","FABP4")
p_vln <- VlnPlot(conbined, features = selected_marker_genes, group.by = "seurat_clusters",
                 pt.size = 0.1) +
  ggtitle("Violin Plot of Selected Marker Genes") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_vln)
p_feature <- FeaturePlot(conbined, features = selected_marker_genes,
                         cols = c("lightgrey", "blue"), ncol = 2) +
  theme(plot.title = element_text(hjust = 0.5))
print(p_feature)
library(Seurat)
library(ggplot2)
selected_marker_genes <- c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")
conbined$group <- factor(conbined$group, levels = c("Healthy", "NAFLD"))
p_feature <- FeaturePlot(conbined, features = selected_marker_genes,
                         cols = c("lightgrey", "blue"), ncol = 2,
                         split.by = "group") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_feature)
selected_marker_genes <- c("IGF1", "MYH11", "HYOU1", "SPATA18", "SCD")
p_feature <- FeaturePlot(
  object = conbined,
  features = selected_marker_genes,
  cols = c("lightgrey", "blue"),
  ncol = 2
) +
  theme(plot.title = element_text(hjust = 0.5))
print(p_feature)
DimPlot(conbined, group.by = "group") + ggtitle("Cell Groups")
print(p_feature)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
selected_genes <- c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")
expression_data <- FetchData(conbined, vars = c(selected_genes, "group"))
long_data <- expression_data %>%
  pivot_longer(cols = all_of(selected_genes),
               names_to = "gene",
               values_to = "expression")
average_expression <- long_data %>%
  group_by(group, gene) %>%
  summarise(mean_expression = mean(expression), .groups = "keep")
p_values_bar <- long_data %>%
  group_by(gene) %>%
  summarise(p_value = wilcox.test(expression ~ group)$p.value)
average_expression <- average_expression %>%
  left_join(p_values_bar, by = "gene")
bar_plot <- ggplot(average_expression, aes(x = gene, y = mean_expression, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(title = "Bar Plot: Mean Gene Expression by Group",
       x = "Gene",
       y = "Mean Expression Level") +
  theme(plot.title = element_text(hjust = 0.5))
for (i in 1:nrow(p_values_bar)) {
  bar_plot <- bar_plot +
    annotate("text",
             x = i,
             y = max(average_expression$mean_expression) * 1.1,
             label = paste("p =", formatC(p_values_bar$p_value[i], format = "e", digits = 2)),
             size = 4, hjust = 0.5)
}
print(bar_plot)
conbined<- JoinLayers(conbined)
macrophages_subset <- subset(conbined, celltype == "Hepatocytes")
Idents(macrophages_subset) <- macrophages_subset$group
markers_diff <- FindMarkers(
  object = macrophages_subset,
  ident.1 = "Healthy",
  ident.2 = "NAFLD",
  features = "NR4A2",
  min.pct = 0.1,
  logfc.threshold = 0
)
print(markers_diff)
p_vln_igf1 <- VlnPlot(macrophages_subset,
                      features = "NR4A2",
                      group.by = "group",
                      pt.size = 0.1) +
  ggtitle("Expression of NR4A2 in Macrophages between Groups") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_vln_igf1)
message("Visualizing IGF1 expression...")
if ("IGF1" %in% rownames(conbined)) {
  p_igf1_umap <- FeaturePlot(conbined, features = "IGF1",
                             pt.size = 0.5, order = TRUE) +
    ggtitle("IGF1 Expression on UMAP")
  p_igf1_violin <- VlnPlot(conbined, features = "IGF1",
                           group.by = "celltype", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("IGF1 Expression by Cell Type")
  p_igf1_violin_group <- VlnPlot(conbined, features = "IGF1",
                                 group.by = "group", pt.size = 0) +
    ggtitle("IGF1 Expression by Group (Healthy vs NAFLD)")
  p_igf1_violin_split <- VlnPlot(conbined, features = "IGF1",
                                 group.by = "celltype",
                                 split.by = "group", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("IGF1 Expression by Cell Type and Group")
  p_igf1_dot <- DotPlot(conbined, features = "IGF1",
                        group.by = "celltype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("IGF1 Expression Level by Cell Type")
  p_igf1_split <- FeaturePlot(conbined, features = "IGF1",
                              split.by = "group", pt.size = 0.5, order = TRUE) +
    ggtitle("IGF1 Expression by Group")
} else {
  message("IGF1 gene not found in the dataset")
  igf_genes <- grep("^IGF", rownames(conbined), value = TRUE)
  if (length(igf_genes) > 0) {
    message("Found these IGF-related genes: ", paste(igf_genes, collapse = ", "))
  }
}
message("Creating basic visualizations...")
p_sample <- DimPlot(conbined, reduction = "umap", group.by = "sample_id", label = TRUE) +
  ggtitle("UMAP by Sample")
p_cluster <- DimPlot(conbined, reduction = "umap", label = TRUE) +
  ggtitle("UMAP by Cluster")
p_celltype <- DimPlot(conbined, reduction = "umap", group.by = "celltype", label = TRUE) +
  ggtitle("UMAP by Cell Type")
p_group <- DimPlot(conbined, reduction = "umap", group.by = "group") +
  ggtitle("UMAP by Group (Healthy vs NAFLD)")
p_split <- DimPlot(conbined, reduction = "umap", group.by = "celltype",
                   split.by = "group", ncol = 2) +
  ggtitle("Cell Types by Group")
message("Saving results...")
save(conbined, file = "conbined_annotated.rdata")
pdf("Complete_Analysis_Results.pdf", width = 12, height = 10)
print(p_sample)
print(p_cluster)
print(p_celltype)
print(p_group)
print(p_split)
print(p_dot)
if ("IGF1" %in% rownames(conbined)) {
  print(p_igf1_umap)
  print(p_igf1_violin)
  print(p_igf1_violin_group)
  print(p_igf1_violin_split)
  print(p_igf1_dot)
  print(p_igf1_split)
}
p_qc1 <- VlnPlot(conbined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "sample_id", pt.size = 0, ncol = 3)
print(p_qc1)
p_qc2 <- VlnPlot(conbined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "group", pt.size = 0, ncol = 3)
print(p_qc2)
dev.off()
message("Analyzing cell proportions...")
library(tidyr)
cell_proportions <- conbined@meta.data %>%
  group_by(sample_id, group, celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_id) %>%
  mutate(proportion = count / sum(count) * 100)
p_prop_sample <- ggplot(cell_proportions,
                        aes(x = sample_id, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportions by Sample",
       x = "Sample", y = "Percentage (%)")
p_prop_group <- ggplot(cell_proportions,
                       aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  labs(title = "Cell Type Proportions by Group",
       x = "Group", y = "Proportion")
p_prop_box <- ggplot(cell_proportions,
                     aes(x = celltype, y = proportion, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportions: Healthy vs NAFLD",
       x = "Cell Type", y = "Percentage (%)")
pdf("Cell_Proportion_Analysis.pdf", width = 12, height = 8)
print(p_prop_sample)
print(p_prop_group)
print(p_prop_box)
dev.off()
message("Analyzing cell proportions...")
library(tidyr)
library(broom)
cell_proportions <- conbined@meta.data %>%
  group_by(sample_id, group, celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_id) %>%
  mutate(proportion = count / sum(count) * 100)
p_values <- cell_proportions %>%
  group_by(celltype) %>%
  summarise(p_value = wilcox.test(proportion ~ group)$p.value)
cell_proportions <- cell_proportions %>%
  left_join(p_values, by = "celltype")
p_prop_sample <- ggplot(cell_proportions,
                        aes(x = sample_id, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportions by Sample",
       x = "Sample", y = "Percentage (%)")
p_prop_group <- ggplot(cell_proportions,
                       aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +
  theme_classic() +
  labs(title = "Cell Type Proportions by Group",
       x = "Group", y = "Proportion")
print(p_prop_group)
p_prop_group <- ggplot(cell_proportions,
                       aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  labs(title = "Cell Type Proportions by Group",
       x = "Group", y = "Proportion")
p_prop_box <- ggplot(cell_proportions,
                     aes(x = celltype, y = proportion, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportions: Healthy vs NAFLD",
       x = "Cell Type", y = "Percentage (%)")
print(p_prop_box)
for (i in unique(cell_proportions$celltype)) {
  p_value <- p_values$p_value[p_values$celltype == i]
  p_prop_box <- p_prop_box +
    annotate("text",
             x = i,
             y = max(cell_proportions$proportion[cell_proportions$celltype == i]) + 1,
             label = paste("p =", formatC(p_value, format = "f", digits = 5)),
             size = 4, hjust = 0.5)
}
pdf("Cell_Proportion_Analysis_with_p_values.pdf", width = 12, height = 8)
print(p_prop_sample)
print(p_prop_group)
print(p_prop_box)
dev.off()
message("Performing differential expression analysis...")
Idents(conbined) <- "celltype"
celltypes <- unique(conbined$celltype)
de_results <- list()
for (ct in celltypes) {
  message("Analyzing: ", ct)
  Idents(conbined) <- "celltype"
  subset_cells <- subset(conbined, idents = ct)
  Idents(subset_cells) <- "group"
  cell_counts <- table(subset_cells$group)
  message("Cell counts: Healthy=", cell_counts["Healthy"], ", NAFLD=", cell_counts["NAFLD"])
  if (cell_counts["Healthy"] >= 3 & cell_counts["NAFLD"] >= 3) {
    de_genes <- FindMarkers(subset_cells,
                            ident.1 = "NAFLD",
                            ident.2 = "Healthy",
                            min.pct = 0.1,
                            logfc.threshold = 0.25)
    de_genes$gene <- rownames(de_genes)
    de_genes$celltype <- ct
    de_results[[ct]] <- de_genes
    if ("IGF1" %in% rownames(de_genes)) {
      message("IGF1 found in DE results for ", ct,
              ": log2FC = ", round(de_genes["IGF1", "avg_log2FC"], 3),
              ", p_val_adj = ", format(de_genes["IGF1", "p_val_adj"], scientific = TRUE))
    }
  }
}
message("Generating analysis report...")
dir.create("analysis_results", showWarnings = FALSE)
save(conbined, de_results, cell_proportions,
     file = "analysis_results/key_analysis_results.rdata")
write.csv(cell_proportions, "analysis_results/cell_proportions.csv", row.names = FALSE)
if (length(de_results) > 0) {
  all_de <- do.call(rbind, de_results)
  write.csv(all_de, "analysis_results/all_differential_expression.csv", row.names = FALSE)
}
cat("=== ANALYSIS COMPLETED SUCCESSFULLY ===\n\n",
    "Total cells:", ncol(conbined), "\n",
    "Cell types:", paste(unique(conbined$celltype), collapse = ", "), "\n",
    "Groups: Healthy (", sum(conbined$group == "Healthy"), " cells), ",
    "NAFLD (", sum(conbined$group == "NAFLD"), " cells)\n",
    "Samples:", paste(unique(conbined$sample_id), collapse = ", "), "\n\n",
    "Key files generated:\n",
    "- conbined_annotated.rdata: Annotated Seurat object\n",
    "- Complete_Analysis_Results.pdf: All visualizations\n",
    "- Cell_Proportion_Analysis.pdf: Cell proportion analysis\n",
    "- analysis_results/: Directory with detailed results\n",
    file = "analysis_summary.txt")
message("Analysis completed!")
message("Check the following files:")
message("1. Complete_Analysis_Results.pdf - All visualizations")
message("2. Cell_Proportion_Analysis.pdf - Cell proportion analysis")
message("3. analysis_results/ - Detailed results directory")
message("4. analysis_summary.txt - Analysis summary")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle3")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(monocle3)
library(SeuratWrappers)
library(monocle3)
set.seed(123)
load("conbined_annotated.rdata")
celltypes_to_analyze <- c("Macrophages" ,"Monocytes")
subset_seurat <- subset(conbined, subset = celltype %in% celltypes_to_analyze)
cat("Selected cell types for trajectory analysis:\n")
print(table(subset_seurat$celltype))
cat("Total cells:", ncol(subset_seurat), "\n")
subset_seurat <- NormalizeData(subset_seurat)
subset_seurat <- FindVariableFeatures(subset_seurat, nfeatures = 2000)
subset_seurat <- ScaleData(subset_seurat)
cds <- as.cell_data_set(subset_seurat)
cds <- cluster_cells(cds)
expression_matrix <- GetAssayData(subset_seurat, assay = "RNA", slot = "counts")
cell_metadata <- subset_seurat@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)
cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds <- cluster_cells(cds)
cds <- learn_graph(cds,verbose = T,learn_graph_control = list(minimal_branch_len=80))
plot1 <- plot_cells(cds,
                    color_cells_by = "celltype",
                    label_groups_by_cluster = FALSE,
                    label_leaves = FALSE,
                    label_branch_points = F,
                    graph_label_size = 3) +
  ggtitle("Trajectory by Cell Type")
plot2 <- plot_cells(cds,
                    color_cells_by = "group",
                    label_groups_by_cluster = FALSE,
                    label_leaves = FALSE,
                    label_branch_points = F,
                    graph_label_size = 3) +
  ggtitle("Trajectory by Group (Healthy vs NAFLD)")
plot3 <- plot_cells(cds,
                    color_cells_by = "seurat_clusters",
                    label_cell_groups = FALSE,
                    label_leaves = TRUE,
                    label_branch_points = F,
                    graph_label_size = 3) +
  ggtitle("Trajectory with Branch Points")
print(plot1)
print(plot2)
print(plot3)
plot_cells(cds,
           color_cells_by = "group",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)
cds <- order_cells(cds)
healthy_cells <- colnames(cds)[cds$group == "Healthy"]
cds <- order_cells(cds, root_cells = healthy_cells)
plot_ordered <- plot_cells(cds,
                           color_cells_by = "pseudotime",
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           graph_label_size = 3) +
  scale_color_viridis_c() +
  ggtitle("Pseudotime Trajectory")
print(plot_ordered)
pseudotime_values <- pseudotime(cds)
subset_seurat$pseudotime <- pseudotime_values
conbined$pseudotime <- NA
conbined$pseudotime[colnames(subset_seurat)] <- pseudotime_values
p_pseudo_umap <- FeaturePlot(conbined,
                             features = "pseudotime",
                             cols = c("lightgrey", "blue"),
                             pt.size = 0.5) +
  ggtitle("Pseudotime Distribution on UMAP") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_pseudo_umap)
deg_cells <- subset_seurat[, !is.na(subset_seurat$pseudotime)]
deg_cells <- NormalizeData(deg_cells)
deg_cells$pseudotime_bin <- cut(deg_cells$pseudotime, breaks = 5)
pseudotime_markers <- FindAllMarkers(deg_cells,
                                     assay = "RNA",
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25,
                                     group.by = "pseudotime_bin")
top_pseudotime_genes <- pseudotime_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
print("Top genes along pseudotime:")
print(top_pseudotime_genes)
p_violin <- plot_genes_violin(cds_subset, group_cells_by="embryo.time.bin", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
print(p_violin)
key_genes <- c("IGF1", "MYH11", "HYOU1", "SPATA18", "SCD")
available_genes <- key_genes[key_genes %in% rownames(deg_cells)]
cat("Available genes for pseudotime expression:", available_genes, "\n")
if (length(available_genes) > 0) {
  pseudotime_expression_plot <- FeaturePlot(deg_cells,
                                            features = available_genes,
                                            cols = c("lightgrey", "red"),
                                            combine = FALSE)
  for (i in 1:length(pseudotime_expression_plot)) {
    pseudotime_expression_plot[[i]] <- pseudotime_expression_plot[[i]] +
      ggtitle(paste(available_genes[i], "Expression along Pseudotime"))
  }
  conbined_expression_plot <- patchwork::wrap_plots(pseudotime_expression_plot, ncol = 2)
  print(conbined_expression_plot)
}
pseudotime_df <- data.frame(
  pseudotime = subset_seurat$pseudotime,
  group = subset_seurat$group,
  celltype = subset_seurat$celltype
)
pseudotime_boxplot <- ggplot(pseudotime_df, aes(x = group, y = pseudotime, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Group",
       x = "Group",
       y = "Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5))
print(pseudotime_boxplot)
healthy_pseudotime <- pseudotime_df$pseudotime[pseudotime_df$group == "Healthy"]
nafld_pseudotime <- pseudotime_df$pseudotime[pseudotime_df$group == "NAFLD"]
if (length(healthy_pseudotime) > 1 & length(nafld_pseudotime) > 1) {
  t_test_result <- t.test(healthy_pseudotime, nafld_pseudotime)
  cat("T-test result for pseudotime difference between groups:\n")
  print(t_test_result)
}
saveRDS(conbined, file = "conbined_with_pseudotime.rds")
saveRDS(cds, file = "monocle3_cds.rds")
cat("Pseudotime analysis completed successfully!\n")
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggplot2)
celltypes_to_analyze <- c("Macrophages")
subset_seurat <- subset(conbined, subset = celltype %in% celltypes_to_analyze)
cat("Selected cell types for trajectory analysis:\n")
print(table(subset_seurat$celltype))
cat("Total cells:", ncol(subset_seurat), "\n")
subset_seurat <- NormalizeData(subset_seurat)
subset_seurat <- FindVariableFeatures(subset_seurat, nfeatures = 2000)
subset_seurat <- ScaleData(subset_seurat)
cds <- as.cell_data_set(subset_seurat)
cds <- cluster_cells(cds)
expression_matrix <- GetAssayData(subset_seurat, assay = "RNA", slot = "counts")
cell_metadata <- subset_seurat@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)
cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds <- cluster_cells(cds)
cds <- learn_graph(cds, verbose = TRUE, learn_graph_control = list(minimal_branch_len = 80))
p1 <- plot_cells(cds, color_cells_by = "group", cell_size = 0.5, group_label_size = 5)
print(p1)
cds <- cluster_cells(cds, resolution = 1e-5)
plot_cells(cds)
plot_cells(cds, color_cells_by = "partition", group_cells_by = "partition")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(conbined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = 'louvain')
cds <- learn_graph(cds, verbose = TRUE, learn_graph_control = list(minimal_branch_len = 80))
p <- plot_cells(cds, color_cells_by = "group",
                label_leaves = FALSE, label_branch_points = FALSE)
print(p)
ggsave("Monocle3_Epi_Cancer.pdf", p, width = 5, height = 4)
plot_cells(cds,
           color_cells_by = "group",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 4)
cds_gene <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
write.csv(cds_gene, file = "cds_gene.csv")
cds_gene <- read.csv("cds_gene.csv", header = TRUE, row.names = 1)
res_ids <- row.names(subset(cds_gene, q_value < 0.01))
gene_module_df <- find_gene_modules(cds[res_ids, ],
                                    resolution = 10^seq(-6, -1))
write.csv(gene_module_df, file = "gene_module_df.csv")
cell_group_df <- tibble::tibble(cell = row.names(colData(cds)),
                                cell_group = colData(cds)$group)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale = "column",
                   clustering_method = "ward.D2")
genes_sig <- cds_gene %>% top_n(n = 4, morans_I) %>% pull(gene_short_name) %>% as.character()
cds <- order_cells(cds)
plot_genes_in_pseudotime(cds[genes_sig, ], color_cells_by = "group", min_expr = 0.5, ncol = 2)
plot_genes_in_pseudotime(cds[genes_sig, ], color_cells_by = "pseudotime", min_expr = 0.5, ncol = 2)
plot_cells(cds, genes = genes_sig, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)
library(viridis)
key_genes <- c("IGF1", "MYH11", "HYOU1", "SPATA18", "SCD")
plot_cells(cds,
           genes = key_genes,
           label_cell_groups = FALSE,
           show_trajectory_graph = TRUE,
           cell_size = 1,
           trajectory_graph_color = "black",
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE) + scale_color_viridis(option = "inferno")
genes_sig <- c("IGF1", "MYH11", "HYOU1", "SPATA18", "SCD")
cds <- order_cells(cds)
plot_genes_in_pseudotime(cds[genes_sig, ], color_cells_by = "group", min_expr = 0.5, ncol = 2)
plot_genes_in_pseudotime(cds[genes_sig, ], color_cells_by = "pseudotime", min_expr = 0.5, ncol = 2)
plot_cells(cds, genes = genes_sig, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)
conbined=conbined
data.input <- GetAssayData(conbined,   slot = "data")
identity <- subset(conbined@meta.data, select = "celltype")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
save(cellchat,file ="cellchat.rdata" )
df.net.1 <- subsetCommunication(cellchat,slot.name = "netP")
df.net.2 <- subsetCommunication(cellchat )
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("IL16", "MK"))
unique(cellchat@netP$pathways)
unique(cellchat@var.features[["features"]])
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
levels(cellchat@idents)
vertex.receiver = c(1, 2)
cellchat@netP$pathways
pathways.show <- "CCL"
?netVisual_aggregate
vertex.receiver = seq(1,2)
netVisual_aggregate(cellchat, signaling = "MK",
                    vertex.receiver = c(1,2,3),layout="hierarchy")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MK", layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MK", vertex.receiver = c(1 ),
                    layout = "chord", vertex.size = groupSize)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "MK", color.heatmap = "Reds")
pathways.show <- "CCL"
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.MK <- extractEnrichedLR(cellchat, signaling = "CCL", geneLR.return = FALSE)
LR.show <- pairLR.MK[2,]
vertex.receiver = seq(1,2)
netVisual_individual(cellchat, signaling ="CCL"  ,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
netVisual_individual(cellchat, signaling ="CCL"  , pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling ="CCL"  , pairLR.use = LR.show, layout = "chord")
pathway.show.all=cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver=c(1,2,3,4)
getwd()
setwd(".")
for (i in 1:length(pathway.show.all)) {
  netVisual(cellchat,signaling = pathway.show.all[i],out.format = c("pdf"),
            vertex.receiver=vertex.receiver,layout="circle")
  plot=netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
  ggsave(filename = paste0(pathway.show.all[i],".contribution.pdf"),
         plot=plot,width=6,height=4,dpi=300,units="in")
}
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use = c(1:5), remove.isolate = FALSE)
cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),
                 signaling =  c("MK","CCL"), remove.isolate = FALSE)
pairLR  <- extractEnrichedLR(cellchat, signaling =c("MK","CCL"), geneLR.return = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),pairLR.use =pairLR , remove.isolate = FALSE)
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:5), lab.cex = 0.5,legend.pos.y = 30)
plotGeneExpression(cellchat, signaling = "MK")
plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE)
plotGeneExpression(cellchat, signaling = "MK",type = "dot")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = "MK", width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat,
                                         signaling = c("MK", "PARs"))
gg1 + gg2
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
table(conbined@meta.data$orig.ident )
sc.sp=SplitObject(conbined,split.by = "orig.ident")
sc.11=conbined[,sample(colnames(sc.sp[["sample2"]]),1000)]
sc.3=conbined[,sample(colnames(sc.sp[["sample21"]]),1000)]
cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")
dir.create("compare")
setwd("compare/")
cellchat=cellchat.sc11
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat
cellchat=cellchat.sc3
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat
cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")
rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )
netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )
table(conbined@active.ident)
s.cell=c( "Macrophage", "Tissue_stem_cells","Monocyte")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]
netVisual_circle(count1,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )
netVisual_circle(count2,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )
library(ggplot2)
library(cowplot)
library(ggrepel)
theme_set(theme_cowplot())
conbined <- conbined
conbined$group <- as.character(conbined$group)
deg <- FindMarkers(conbined, ident.1 = "Healthy", ident.2 = "NAFLD",
                   group.by = "group", logfc.threshold = 0.01)
exp <- AverageExpression(conbined, group.by = "group", verbose = FALSE)
exp <- as.data.frame(exp$RNA)
exp <- log1p(exp)
avg.expression <- exp
avg.expression$gene <- rownames(avg.expression)
genes.to.label <- sample(rownames(deg), 20)
p1 <- ggplot(avg.expression, aes(x = Healthy, y = NAFLD)) +
  geom_point() +
  ggtitle("Healthy vs NAFLD") +
  geom_label_repel(data = avg.expression[rownames(avg.expression) %in% genes.to.label,],
                   aes(label = gene),
                   size = 3,
                   show.legend = FALSE)
print(p1)
deg$symbol <- rownames(deg)
logFC_t <- 0
P.Value_t <- 1e-28
deg$change <- ifelse(deg$p_val_adj < P.Value_t & deg$avg_log2FC < 0, "down",
                     ifelse(deg$p_val_adj < P.Value_t & deg$avg_log2FC > 0, "up", "stable"))
library(ggplot2)
library(ggrepel)
p_volcano <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), size = pct.1)) +
  geom_point(alpha = 0.6, aes(color = change), stroke = 1.2) +
  ylab("-log10(Pvalue)") +
  xlab("Log2 Fold Change") +
  scale_color_manual(values = c("green", "grey", "red")) +
  geom_hline(yintercept = -log10(P.Value_t), lty = 2, col = "black", lwd = 0.8) +
  theme_minimal(base_size = 15) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  ggtitle("Volcano Plot: Healthy vs NAFLD") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
data_selected <- deg[deg$symbol == "IGF1", ]
p_volcano + geom_label_repel(data = data_selected, aes(label = symbol), size = 5, color = "black")
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(data.table)
library(tidyverse)
conbined <- readRDS("conbined.rds")
table(conbined@meta.data$celltype)
table(conbined@meta.data$orig.ident)
DimPlot(conbined, group.by = "orig.ident")
DimPlot(conbined, group.by = "celltype")
type <- unique(conbined@meta.data$celltype)
dir.create("deg/")
setwd("deg/")
r.deg <- data.frame()
table(conbined@meta.data$orig.ident)
for (i in 1:length(type)) {
  Idents(conbined) <- "celltype"
  deg <- FindMarkers(conbined, ident.1 = "NAFLD", ident.2 = "Healthy",
                     group.by = "group", subset.ident = type[i])
  write.csv(deg, file = paste0(type[i], 'deg.csv'))
  deg$gene <- rownames(deg)
  deg$celltype <- type[i]
  deg$unm <- i-1
  r.deg <- rbind(deg, r.deg)
}
r.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0, 'Up', 'Down'))
dim(r.deg)
r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.01, 'Highly', 'Lowly'))
r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
r.deg$unm %<>% as.vector(.) %>% as.numeric(.)
celltype_mapping <- r.deg %>%
  dplyr::select(unm, celltype) %>%
  distinct() %>%
  arrange(unm)
highlight_genes <- c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")
r.deg$is_highlight <- ifelse(r.deg$gene %in% highlight_genes, "Highlight", "Normal")
r.deg$thr_signi_highlight <- ifelse(r.deg$is_highlight == "Highlight",
                                    paste0(r.deg$threshold, "_Highlight"),
                                    r.deg$thr_signi)
top_up_label <- r.deg %>%
  subset(., threshold %in% "Up") %>%
  group_by(unm) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  as.data.frame()
top_down_label <- r.deg %>%
  subset(., threshold %in% "Down") %>%
  group_by(unm) %>%
  top_n(n = -5, wt = avg_log2FC) %>%
  as.data.frame()
top_label <- rbind(top_up_label, top_down_label)
background_position <- r.deg %>%
  dplyr::group_by(unm) %>%
  dplyr::summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  as.data.frame()
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black",
                    "Up_Highlight" = "red",
                    "Down_Highlight" = "red")
celltype_count <- length(unique(r.deg$celltype))
cols_cluster <- brewer.pal(n = max(8, celltype_count), name = "Set3")[1:celltype_count]
names(cols_cluster) <- as.character(0:(celltype_count-1))
p <- ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min, ymax = Max),
            fill = "#525252", alpha = 0.1) +
  geom_jitter(data = r.deg, aes(x = unm, y = avg_log2FC, colour = thr_signi_highlight),
              size = 1, position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi,
                     breaks = c("Up_Highly", "Down_Highly", "Up_Lowly", "Down_Lowly"),
                     labels = c("Up Highly", "Down Highly", "Up Lowly", "Down Lowly")) +
  scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                     breaks = seq(0, max(r.deg$unm), 1),
                     labels = celltype_mapping$celltype) +
  geom_text_repel(data = top_label, aes(x = unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                             ymax = 0.4, fill = as.factor(unm)),
            color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "Cell Type", y = "average log2FC") +
  theme_bw()
plot1 <- p + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.text.y = element_text(colour = 'black', size = 14),
  axis.text.x = element_text(colour = 'black', size = 12, angle = 45, hjust = 1, vjust = 1),
  panel.border = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.y = element_line(colour = "black")
)
plot1
ggsave(filename = "deg_pointplot.pdf", plot = plot1, width = 9, height = 6)
highlight_genes_in_data <- r.deg[r.deg$gene %in% highlight_genes, ]
if (nrow(highlight_genes_in_data) > 0) {
  cat("以下突出显示基因在差异分析结果中:\n")
  print(highlight_genes_in_data[, c("gene", "celltype", "avg_log2FC", "p_val_adj")])
} else {
  cat("没有找到突出显示的基因在差异分析结果中\n")
}
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
readRDS("conbined.rds")
conbined=.Last.value
table(conbined@meta.data$celltype)
table(conbined@meta.data$orig.ident)
DimPlot(conbined,group.by = "orig.ident")
DimPlot(conbined,group.by = "celltype")
table(conbined@meta.data$celltype)
table(conbined@meta.data$orig.ident)
type=unique(conbined@meta.data$celltype)
dir.create("deg/")
setwd("deg/")
r.deg=data.frame()
table(conbined@meta.data$orig.ident)
for (i in 1:length(type)) {
  Idents(conbined)="celltype"
  deg=FindMarkers(conbined,ident.1 = "NAFLD",ident.2 = "Healthy",
                  group.by = "group",subset.ident =type[i]   )
  write.csv(deg,file = paste0( type[i],'deg.csv') )
  deg$gene=rownames(deg)
  deg$celltype=type[i]
  deg$unm=i-1
  r.deg=rbind(deg,r.deg)
}
Idents(conbined)="celltype"
deg=FindMarkers(conbined,ident.1 = "NAFLD",ident.2 = "Healthy",
                group.by = "group",subset.ident =type[2]   )
table(r.deg$unm)
r.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0 , 'Up', 'Down'))
dim(r.deg)
r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
r.deg$unm %<>% as.vector(.) %>% as.numeric(.)
top_up_label <- r.deg %>%
  subset(., threshold%in%"Up") %>%
  group_by(unm) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  as.data.frame()
top_down_label <- r.deg %>%
  subset(., threshold %in% "Down") %>%
  group_by(unm) %>%
  top_n(n = -5, wt = avg_log2FC) %>%
  as.data.frame()
top_label <- rbind(top_up_label,top_down_label)
top_label$thr_signi %<>%
  factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))
colnames(r.deg)
background_position <- r.deg %>%
  dplyr::group_by(unm) %>%
  dplyr::summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  as.data.frame()
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cluster_bar_position$unm %<>%
  factor(., levels = c(0:max(as.vector(.))))
cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black")
cols_cluster <- c("0" = "#35978f",
                  "1" = "#8dd3c7",
                  "2" = "#ffffb3",
                  "3" = "#bebada",
                  "4" = "#fb8072",
                  "5" = "#80b1d3",
                  "6" = "#fdb462","7" = "#925bea","8" = "#db5e92"
)
p= ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "#525252", alpha = 0.1) +
  geom_jitter(data = r.deg, aes(x =unm, y = avg_log2FC, colour = thr_signi),
              size = 1,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi) +
  scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                     breaks = seq(0, max(r.deg$unm), 1),
                     label = seq(0, max(r.deg$unm),1)) +
  geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                             ymax = 0.4, fill = unm), color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "Cluster", y = "average log2FC") +
  theme_bw()
plot1 <- p + theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour =c("#89288F","#89288F","#F47D2B","#F47D2B","#FF00FF",'black','black'), size = 14, vjust = 58),
                   panel.border = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.line.y = element_line(colour = "black"))
plot1
ggsave(filename = "deg_pointplot.pdf", plot = plot1, width = 9, height = 6)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
setwd(".")
load("conbined_annotated.rdata")
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)
data.input <- GetAssayData(conbined,   slot = "data")
identity <- subset(conbined@meta.data, select = "celltype")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
target_celltype <- "CellTypeA"
df.communication <- subsetCommunication(cellchat, sources.use = "Hepatocytes")
head(df.communication)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
source_cells <- which(cellchat@idents == target_celltype)
target_cells <- which(cellchat@idents != target_celltype)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 sources.use = source_cells, targets.use = target_cells,
                 weight.scale = T, label.edge = F,
                 title.name = paste("Interactions for", target_celltype))
head(cellchat@net$count)
rownames(cellchat@net$count) <- names(cellchat@net$count)
head(cellchat@net$count)
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
levels(cellchat@idents)
vertex.receiver = c(6)
cellchat@netP$pathways
pathways.show <- "VISFATIN"
?netVisual_aggregate
vertex.receiver = seq(6,8)
netVisual_aggregate(cellchat, signaling = "MIF",
                    vertex.receiver = vertex.receiver,layout="hierarchy")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MIF", layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MIF", vertex.receiver = c(1 ),
                    layout = "chord", vertex.size = groupSize)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "MIF", color.heatmap = "Reds")
pathways.show <- "CCL"
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.MK <- extractEnrichedLR(cellchat, signaling = "CCL", geneLR.return = FALSE)
LR.show <- pairLR.MK[2,]
vertex.receiver = seq(1,2)
netVisual_individual(cellchat, signaling ="CCL"  ,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
netVisual_individual(cellchat, signaling ="CCL"  , pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling ="CCL"  , pairLR.use = LR.show, layout = "chord")
pathway.show.all=cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver=c(6,8)
getwd()
setwd(".")
for (i in 1:length(pathway.show.all)) {
  netVisual(cellchat,signaling = pathway.show.all[i],out.format = c("pdf"),
            vertex.receiver=vertex.receiver,layout="circle")
  plot=netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
  ggsave(filename = paste0(pathway.show.all[i],".contribution.pdf"),
         plot=plot,width=6,height=4,dpi=300,units="in")
}
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:5), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(6,8), targets.use = c(1:5), remove.isolate = FALSE)
cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use =c(1,5), targets.use =c(6,8),
                 signaling =  c("MIF","VISFATIN"), remove.isolate = FALSE)
pairLR  <- extractEnrichedLR(cellchat, signaling =c("MIF","CCL"), geneLR.return = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),pairLR.use =pairLR , remove.isolate = FALSE)
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:5), lab.cex = 0.5,legend.pos.y = 30)
plotGeneExpression(cellchat, signaling = "MK")
plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE)
plotGeneExpression(cellchat, signaling = "MK",type = "dot")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = "MK", width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat,
                                         signaling = c("MK", "PARs"))
gg1 + gg2
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
table(conbined@meta.data$orig.ident )
sc.sp=SplitObject(conbined,split.by = "orig.ident")
sc.11=conbined[,sample(colnames(sc.sp[["sample2"]]),1000)]
sc.3=conbined[,sample(colnames(sc.sp[["sample21"]]),1000)]
cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")
dir.create("compare")
setwd("compare/")
cellchat=cellchat.sc11
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat
cellchat=cellchat.sc3
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat
cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")
rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )
netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )
table(conbined@active.ident)
s.cell=c( "Macrophage", "Tissue_stem_cells","Monocyte")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]
netVisual_circle(count1,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )
netVisual_circle(count2,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )
setwd(".")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
setwd(".")
load("combined_annotated.rdata")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
readRDS("conbined.rds")
combined=.Last.value
table(combined@meta.data$celltype)
table(combined@meta.data$seurat_clusters)
Idents(combined)="celltype"
scRNA.Osteoclastic=subset(combined,ident=c("T Cells"))
data <- GetAssayData(scRNA.Osteoclastic,slot = "count")
data[1:20,1:20]
pd <- new('AnnotatedDataFrame', data = scRNA.Osteoclastic@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "celltype")
plot_cell_trajectory(HSMM, color_by = "group")
plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "celltype")
plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~State, nrow = 3)
plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~orig.ident )
blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")))
plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)
HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by =  "celltype")
genes <- c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")
p1 <- plot_genes_jitter(HSMM[genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(HSMM[genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(HSMM[genes,], color_by = "State")
plotc <- p1|p2|p3
plotc
genes <- c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")
p1 <- plot_genes_jitter(HSMM[genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(HSMM[genes,], grouping = "State", color_by = "State")
plotc <- p1|p2
plotc
to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4")))
cds_subset <- HSMM[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by ="State")
dev.off()
marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("BCL6", "IL10", "SHBG", "NR4A2", "CASP1", "FABP4","MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))
diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 1))
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = T)
plot_cell_trajectory(HSMM, color_by = "State")
BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 7)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
genes <- row.names(subset(fData(HSMM),
                          gene_short_name %in% c( "MEF2C", "CCNB2", "TNNT1")))
plot_genes_branched_pseudotime(HSMM[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)
