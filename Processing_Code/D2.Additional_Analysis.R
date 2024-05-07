library("Seurat")
library("SingleR")
library("dplyr")
library("ggplot2")
library("purrr")
library("gridExtra")
# 0. Setting up enviornment
proj_dir = "/proj/hyejunglab/singlecell/Sangmi"
setwd(proj_dir)
output_dir = "/proj/hyejunglab/singlecell/Sangmi/derived_data"
Process_log = "/proj/hyejunglab/singlecell/Sangmi/Processing_Code/process_log.txt"
seurat_data = readRDS(paste0(output_dir,"/integrated_All.RDS"))
sc1_obj = readRDS(paste0(output_dir,"/sc001_seurat.rds"))
sc2_obj = readRDS(paste0(output_dir,"/sc002_seurat.rds"))

target = sc1_obj
gene_vector = data.frame(cbind(num = 1:dim(target@assays$RNA@features@.Data)[1],gene = rownames(target@assays$RNA@features@.Data)))
marker_gene = strsplit("Gnb2l1, DCX, NKX2-1, GAD1, SOX6, MAF, MAFB, DLX5, DLX2, MEF2C, MAP2, 
LHX6, ERBB4, SST, MKI67, LHX8, SP8, PAX6, TH, TPH1, CHAT, GFAP, MBP, POU5F1, CDKN1C, PEG10, 
PKIA, ENC1, NRP2, NOVA1, ZFHX4, NTM, TSHZ2, PBX3, CALY, C1orf61, DLX6-AS1, ZEB2, TCF4, 
ARX, GAD2, SLC32A1, CXCR4, CCK, VIP, NPY, RELN, CALB1 ,
GBX2, DNM3, NDN",", ")[[1]]

marker_index = gene_vector %>%
  filter(gene %in% marker_gene)
gene_list = marker_index$gene
gene_list[19] = "DLX6_AS1"
gene_list[38] = "NKX2_1"

umap_total = cbind(seurat_data@meta.data$Batch_Name,seurat_data@reductions[["umap.rpca"]]@cell.embeddings)
umap_total = cbind(rownames(umap_total),umap_total)
colnames(umap_total) = c("cell_ID","Batch_Name","umap1","umap2")

mat = target@assays$RNA@layers$dataå
gene_mat = t(data.frame(mat[as.numeric(marker_index[,1]),]))
gene_mat = data.frame(cbind(paste0("SC001_",rownames(target@assays$RNA@cells@.Data)),gene_mat))
colnames(gene_mat) = c("cell_ID",gene_list)
sc1_gene = gene_mat
sc1_mat = merge(umap_total,sc1_gene,by="cell_ID") %>%
  mutate(across(
    .cols = -c("cell_ID","Batch_Name"),
    .fns = as.numeric
  ))


target = sc2_obj
gene_vector = data.frame(cbind(num = 1:dim(target@assays$RNA@features@.Data)[1],gene = rownames(target@assays$RNA@features@.Data)))
marker_index = gene_vector %>%
  filter(gene %in% marker_gene)
mat = target@assays$RNA@layers$data
gene_mat = t(data.frame(mat[as.numeric(marker_index[,1]),]))
gene_mat = data.frame(cbind(paste0("SC002_",rownames(target@assays$RNA@cells@.Data)),gene_mat))
colnames(gene_mat) = c("cell_ID",gene_list)
sc2_gene = gene_mat
sc2_mat = merge(umap_total,sc2_gene,by="cell_ID") %>%
  mutate(across(
    .cols = -c("cell_ID","Batch_Name"),
    .fns = as.numeric
  ))

target = sc1_mat
name = "SC_001"
for (i in 1:length(gene_list)){
  target %>%
    ggplot(aes_string(x="umap1",y="umap2",color=gene_list[i])) +
    geom_point(alpha=0.3,size=0.7) +
    ggtitle(paste(name," UMAP for the gene : ",gene_list[i])) +
    theme(plot.title = element_text(size = 30)) + 
    scale_colour_gradient(low = "#56B1F7",high = "#132B43")
  ggsave(paste0("Result/Gene_plot/",name,"_Gene_plot_",gene_list[i],".png"))
}

figure_list = list()
for (i in 1:length(gene_list)) {
  p <- target %>%
    ggplot(aes_string(x="umap1", y="umap2", color=gene_list[i])) +
    geom_point(alpha=0.3,size=0.5) +
    ggtitle(paste(name," UMAP for the gene:", gene_list[i])) +
    theme(plot.title = element_text(size = 15)) + 
    scale_colour_gradient(low = "#56B1F7", high = "#132B43")
  figure_list[[i]] <- p
}

for (j in seq(1, length(figure_list), by = 6)) {
  plots_to_print <- figure_list[j:min(j+5, length(figure_list))]
  grid_plot <- do.call(grid.arrange, c(plots_to_print, ncol = 3, nrow = 2))
  ggsave(paste0("Result/Gene_grid/",name,"_Gene_grid_", ceiling(j / 6), ".png"), grid_plot, width = 20, height = 10)
}

target = sc2_mat
name = "SC_002"
for (i in 1:length(gene_list)){
  target %>%
    ggplot(aes_string(x="umap1",y="umap2",color=gene_list[i])) +
    geom_point(alpha=0.3,size=0.7) +
    ggtitle(paste(name," UMAP for the gene : ",gene_list[i])) +
    theme(plot.title = element_text(size = 30)) + 
    scale_colour_gradient(low = "#56B1F7",high = "#132B43")
  ggsave(paste0("Result/Gene_plot/",name,"_Gene_plot_",gene_list[i],".png"))
}

figure_list = list()
for (i in 1:length(gene_list)) {
  p <- target %>%
    ggplot(aes_string(x="umap1", y="umap2", color=gene_list[i])) +
    geom_point(alpha=0.3,size=0.5) +
    ggtitle(paste(name," UMAP for the gene:", gene_list[i])) +
    theme(plot.title = element_text(size = 15)) + 
    scale_colour_gradient(low = "#56B1F7", high = "#132B43")
  figure_list[[i]] <- p
}

for (j in seq(1, length(figure_list), by = 6)) {
  plots_to_print <- figure_list[j:min(j+5, length(figure_list))]
  grid_plot <- do.call(grid.arrange, c(plots_to_print, ncol = 3, nrow = 2))
  ggsave(paste0("Result/Gene_grid/",name,"_Gene_grid_", ceiling(j / 6), ".png"), grid_plot, width = 20, height = 10)
}

library(reshape2)

n_gene = length(gene_list)
n_plot = n_gene%/%10
i=0
colnames(sc1_mat[,(5+i*10):(4+(i+1)*10)])
i=1
colnames(sc1_mat[,(5+i*10):(4+(i+1)*10)])
i=2
colnames(sc1_mat[,(5+i*10):(4+(i+1)*10)])
i=3
colnames(sc1_mat[,(5+i*10):(4+(i+1)*10)])
colnames(sc1_mat)

sc1_melt = melt(sc1_mat)

for (i in 0:(n_plot-1)){
  first = 5+i*10 
  if (5+(i+1)*10 < length(colnames(sc1_mat))){
    last = (5+(i+1)*10)
  } else {
    last = length(colnames(sc1_mat))
  }
  sc1_melt = melt(sc1_mat[first:last])
  sc2_melt = melt(sc2_mat[first:last])

  p1 = sc1_melt %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 1, fill = "skyblue", colour = "darkblue") +
    labs(title = "Boxplot of Gene expression", x = "Genes", y = "Expression") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5), # 제목을 중앙 정렬
      axis.text.x = element_text(angle = 45, hjust = 1) # X축 라벨을 기울여서 표시
    ) + 
    ylim(c(0,8))
  p2 = sc2_melt %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 1, fill = "skyblue", colour = "darkblue") +
    labs(title = "Boxplot of Gene expression", x = "Genes", y = "Expression") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5), # 제목을 중앙 정렬
      axis.text.x = element_text(angle = 45, hjust = 1) # X축 라벨을 기울여서 표시
    )+ 
    ylim(c(0,8))
  ggsave(paste0("Result/Gene_Boxplot/Gene_Boxplot_",i,".png"), grid.arrange(p1,p2,nrow=2), width = 20, height = 10)

  p1 = sc1_melt %>%
    filter(value != 0) %>%
    ggplot(aes(x= variable,y=value)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 1, fill = "skyblue", colour = "darkblue") +
    labs(title = "Boxplot of Gene expression (Non Zero)", x = "Genes", y = "Expression") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5), # 제목을 중앙 정렬
      axis.text.x = element_text(angle = 45, hjust = 1) # X축 라벨을 기울여서 표시
    )+ 
    ylim(c(0,8))
  p2 = sc2_melt %>%
    filter(value != 0) %>%
    ggplot(aes(x= variable,y=value)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 1, fill = "skyblue", colour = "darkblue") +
    labs(title = "Boxplot of Gene expression (Non Zero)", x = "Genes", y = "Expression") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5), # 제목을 중앙 정렬
      axis.text.x = element_text(angle = 45, hjust = 1) # X축 라벨을 기울여서 표시
    )+ 
    ylim(c(0,8))
  ggsave(paste0("Result/Gene_Boxplot/Gene_nonzero_Boxplot_",i,".png"), grid.arrange(p1,p2,nrow=2), width = 20, height = 10)
}
  

