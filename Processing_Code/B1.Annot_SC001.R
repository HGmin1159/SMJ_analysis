library("Seurat")
library("SingleR")
library("scater")

# 0. Setting up enviornment
proj_dir = "/proj/hyejunglab/singlecell/Sangmi"
setwd(proj_dir)
output_dir = "/proj/hyejunglab/singlecell/Sangmi/derived_data"
Process_log = "/proj/hyejunglab/singlecell/Sangmi/Processing_Code/process_log.txt"

cat("Done : Setting up enviornment \n",file = Process_log)

# 1. Labeling and DR the data

ref_seurat_path = paste0(proj_dir,"/Data/GSE208672_Seurat_allsamples.rds")
ref_seurat_obj = readRDS(ref_seurat_path)

sc001_cellranger_path = "/proj/hyejunglab/singlecell/Sangmi/Data/SC001_cellranger_count_outs/filtered_feature_bc_matrix"
sc001_seurat_obj = Read10X(sc001_cellranger_path)
sc001_seurat_obj = CreateSeuratObject(counts = sc001_seurat_obj) 

ref.input <- as.SingleCellExperiment(ref_seurat_obj)
ref.input <- logNormCounts(ref.input)

DefaultAssay(sc001_seurat_obj) <- "RNA"
sc001_seurat_obj <- JoinLayers(sc001_seurat_obj)
test.input <- Seurat::as.SingleCellExperiment(sc001_seurat_obj)

pred <- SingleR(test = test.input, ref = ref.input, labels = ref.input$celltype, de.method = "wilcox", assay.type.test = "counts", assay.type.ref = "logcounts")

saveRDS(pred, file = paste0(output_dir,"/sc001_labeled.rds"))

sc001_seurat_obj[["SingleR_anno_subclass"]] <- pred$labels

sc001_seurat_obj <- NormalizeData(sc001_seurat_obj)
sc001_seurat_obj <- FindVariableFeatures(sc001_seurat_obj)
sc001_seurat_obj <- ScaleData(sc001_seurat_obj)
sc001_seurat_obj <- RunPCA(sc001_seurat_obj)

saveRDS(sc001_seurat_obj, file = paste0(output_dir,"/sc001_seurat.rds"))


cat("Done : Labeling and DR \n",file = Process_log,append = TRUE)
