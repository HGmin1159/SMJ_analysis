library("Seurat")
library("SingleR")
library("scater")

# 0. Setting up enviornment
proj_dir = "/proj/hyejunglab/singlecell/Sangmi"
setwd(proj_dir)
output_dir = "/proj/hyejunglab/singlecell/Sangmi/derived_data"
Process_log = "/proj/hyejunglab/singlecell/Sangmi/Processing_Code/process_log.txt"


# 2. Integrating Ref and Test
ref_seurat_path = paste0(proj_dir,"/Data/GSE208672_Seurat_allsamples.rds")
ref_seurat_obj = readRDS(ref_seurat_path)

seurat_obj2 = readRDS(paste0(output_dir,"/sc002_seurat.rds"))
ref_seurat_obj[["Batch_Name"]] = "Reference"
ref_seurat_obj@meta.data[["orig.ident"]] = "Ref"
seurat_obj2[["Batch_Name"]] = "SC002"
seurat_obj2@meta.data[["orig.ident"]] = "SC002"

obj = merge(ref_seurat_obj,seurat_obj2)
obj = NormalizeData(obj)
obj = FindVariableFeatures(obj)
obj = ScaleData(obj)
obj = RunPCA(obj)

obj <- IntegrateLayers(object = obj, method = CCAIntegration,
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "rpca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

saveRDS(obj,paste0(output_dir,"/integrated_SC002.RDS"))

cat("Done : Integrating Ref and Test \n",file = Process_log,append = TRUE)
