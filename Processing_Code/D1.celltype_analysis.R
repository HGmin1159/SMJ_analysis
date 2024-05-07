library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
working_dir = "/proj/hyejunglab/singlecell/Sangmi"
setwd(working_dir)
# seurat_obj = readRDS(paste0(working_dir,"/derived_data/integrated_All.RDS"))

# meta = seurat_obj@meta.data
# head(meta)
# rpca_embeddings <- seurat_obj@reductions[["integrated.rpca"]]@cell.embeddings
# head(rpca_embeddings)
# int_data <- merge(meta, rpca_embeddings, by = "row.names", all.x = TRUE)
# umap_embeddings <- seurat_obj@reductions[["umap.rpca"]]@cell.embeddings
# head(umap_embeddings)
# head(int_data)
# row.names(int_data) = int_data$Row.names
# int_data = int_data[,-1]

# int_data <- merge(int_data, umap_embeddings, by = "row.names", all.x = TRUE)
# head(int_data)
# colnames(int_data)

# meta_data_colnames <- colnames(seurat_obj@meta.data)
# rpca_colnames <- paste0("PCA_", 1:ncol(rpca_embeddings))
# umap_colnames <- paste0("UMAP_", 1:ncol(umap_embeddings))
# colnames(int_data) <- c("Cell_ID",meta_data_colnames, rpca_colnames, umap_colnames)
# int_data[int_data$orig.ident == "SC001","celltype"] = int_data[int_data$orig.ident == "SC001","SingleR_anno_subclass"]
# int_data[int_data$orig.ident == "SC002","celltype"] = int_data[int_data$orig.ident == "SC002","SingleR_anno_subclass"]
# write.csv(int_data, "/proj/hyejunglab/singlecell/Sangmi/derived_data/DR_embedded_table.csv")

data = read.csv("/proj/hyejunglab/singlecell/Sangmi/derived_data/DR_embedded_table.csv")
colnames(data)

sc1_Data = data[data$orig.ident == "SC001",]
sc2_Data = data[data$orig.ident == "SC002",]
ref_Data = data[data$orig.ident == "Ref", ]
data[data$orig.ident=="SC001",]$process = "SC001"
data[data$orig.ident=="SC002",]$process = "SC002"


table(data$orig.ident,useNA="always")
table(data$Batch_Name,useNA="always")
table(data$process,useNA="always")
table(ref_Data["celltype"])

table(sc1_Data["celltype"])
table(sc2_Data["celltype"])


x_range <- range(data$UMAP_1, na.rm = TRUE)
y_range <- range(data$UMAP_2, na.rm = TRUE)

png("Result/umap_celltype.png", width = 1800, height = 600)
p1 = data %>%
    ggplot(aes(x=UMAP_1,y=UMAP_2,group=celltype,color=celltype)) +
    geom_point(size=0.5,alpha = 0.2)+
    ggtitle("Celltype Distribution with UMAP \n for Ref + SC001 + SC002")+
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
p2 = sc1_Data %>%
    ggplot(aes(x=UMAP_1,y=UMAP_2,group=celltype,color=celltype)) +
    geom_point(size=0.5,alpha = 0.2)+
    ggtitle("Celltype Distribution with UMAP \n for SC001")+
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
p3 = sc2_Data %>%
    ggplot(aes(x=UMAP_1,y=UMAP_2,group=celltype,color=celltype)) +
    geom_point(size=0.5,alpha = 0.2)+
    ggtitle("Celltype Distribution with UMAP \n for SC002")+
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

colors <- c("SC001"= "red","SC002" = "blue","Day14"  = "green", "Day0" = "orange",  "Sorted" = "purple", "Unsorted"  = "brown")

png("Result/umap_Process.png", width = 1800, height = 600)
p1 = data %>% 
    ggplot(aes(x=UMAP_1, y=UMAP_2, group=process, color=process)) +
    geom_point(size=0.5,alpha = 0.2) +
    ggtitle("Process Distribution with UMAP \n for Ref+SC001+SC002") +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
p2 = data %>% 
    filter(process %in% c("SC001", "SC002", "Day14")) %>%
    ggplot(aes(x=UMAP_1, y=UMAP_2, group=process, color=process)) +
    geom_point(size=0.5,alpha = 0.2) +
    ggtitle("Process Distribution with UMAP \n for Day14+SC001+SC002") +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
p3 = data %>% 
    filter(process %in% c("Day0", "Day14", "Sorted", "Unsorted")) %>%
    ggplot(aes(x=UMAP_1, y=UMAP_2, group=process, color=process)) +
    geom_point(size=0.5,alpha = 0.2) +
    ggtitle("Process Distribution with UMAP \n for Ref") +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
grid.arrange(p1, p2, p3, nrow=1)
dev.off()

x_range1 = c(-2.5,10)
y_range1 = c(-5,15)

png("Result/umap_Process_zoom.png", width = 900, height = 1000)
data %>% 
    filter(process %in% c("SC001", "SC002", "Day14")) %>%
    ggplot(aes(x=UMAP_1, y=UMAP_2, group=process, color=process)) +
    geom_point(size=1,alpha = 0.2) +
    ggtitle("Process Distribution with UMAP \n for Day14+SC001+SC002") +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range1) +
    ylim(y_range1) +
    scale_color_manual(values = colors)
dev.off()


manova_data = rbind(sc1_Data,sc2_Data)
group_A_pca1 <- manova_data %>% filter(Batch_Name == "SC001") %>% pull(PCA_1)
group_B_pca1 <- manova_data %>% filter(Batch_Name == "SC002") %>% pull(PCA_1)
group_A_pca2 <- manova_data %>% filter(Batch_Name == "SC001") %>% pull(PCA_2)
group_B_pca2 <- manova_data %>% filter(Batch_Name == "SC002") %>% pull(PCA_2)
group_A_umap1 <- manova_data %>% filter(Batch_Name == "SC001") %>% pull(UMAP_1)
group_B_umap1 <- manova_data %>% filter(Batch_Name == "SC002") %>% pull(UMAP_1)
group_A_umap2 <- manova_data %>% filter(Batch_Name == "SC001") %>% pull(UMAP_2)
group_B_umap2 <- manova_data %>% filter(Batch_Name == "SC002") %>% pull(UMAP_2)



png("Result/umap_Batch.png", width = 1200, height = 600)
p1 = data %>% 
    filter(process %in% c("SC001", "SC002")) %>%
    ggplot(aes(x=PCA_1, y=PCA_2, group=process, color=process)) +
    geom_point(size=0.5,alpha = 0.2) +
    ggtitle("Batch Distribution with PCA") +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
p2 = data %>% 
    filter(process %in% c("SC001", "SC002")) %>%
    ggplot(aes(x=UMAP_1, y=UMAP_2, group=process, color=process)) +
    geom_point(size=0.5,alpha = 0.2) +
    ggtitle("Batch Distribution with UMAP") +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    xlim(x_range) +
    ylim(y_range) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))
grid.arrange(p1,p2,nrow=1)
dev.off()


png("Result/density_plot.png", width = 1500, height = 1500)
p1 = ggplot() +
    geom_density(data = data.frame(value = group_A_pca1, Batch_Name = "SC001"), aes(x = value, fill = "SC001"), alpha = 0.5) +
    geom_density(data = data.frame(value = group_B_pca1, Batch_Name = "SC002"), aes(x = value, fill = "SC002"), alpha = 0.5) +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    ggtitle("Density Plot for PCA1")
p2 = ggplot() +
    geom_density(data = data.frame(value = group_A_pca2, Batch_Name = "SC001"), aes(x = value, fill = "SC001"), alpha = 0.5) +
    geom_density(data = data.frame(value = group_B_pca2, Batch_Name = "SC002"), aes(x = value, fill = "SC002"), alpha = 0.5) +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    ggtitle("Density Plot for PCA2")
p3 = ggplot() +
    geom_density(data = data.frame(value = group_A_umap1, Batch_Name = "SC001"), aes(x = value, fill = "SC001"), alpha = 0.5) +
    geom_density(data = data.frame(value = group_B_umap1, Batch_Name = "SC002"), aes(x = value, fill = "SC002"), alpha = 0.5) +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    ggtitle("Density Plot for UMAP1")
p4 = ggplot() +
    geom_density(data = data.frame(value = group_A_umap2, Batch_Name = "SC001"), aes(x = value, fill = "SC001"), alpha = 0.5) +
    geom_density(data = data.frame(value = group_B_umap2, Batch_Name = "SC002"), aes(x = value, fill = "SC002"), alpha = 0.5) +
    theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15))  +
    ggtitle("Density Plot for UMAP2")
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()
