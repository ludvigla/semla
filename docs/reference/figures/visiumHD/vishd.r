library(semla)
library(tidyverse)

st.dir <- "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/binned_outputs"
res <- paste(0, c(8, 16), "um", sep = "")
all.res <- paste(res, collapse = "|")
# retrieve the directories for each resolution
res.dir <- list.dirs(st.dir, recursive = FALSE) %>%
    str_subset(., pattern = all.res)
# build the infoTable for all resolutions
infoTable <- data.frame(
    samples = list.files(res.dir,
        full.names = TRUE, recursive = TRUE,
        pattern = paste0("^filtered_feature.+.h5$")
    ),
    spotfiles = list.files(res.dir,
        full.names = TRUE, recursive = TRUE,
        pattern = "parquet$|positions.csv$"
    ),
    imgs = list.files(res.dir,
        recursive = TRUE,
        full.names = TRUE, pattern = "hires"
    ) %>%
        str_subset(all.res),
    json = list.files(st.dir,
        recursive = TRUE,
        full.names = TRUE, pattern = "^scalefactors"
    ) %>%
        str_subset(all.res),
    resolution = res,
    sample_ID = "mousebrain"
)

se.hd <- ReadVisiumData(infoTable)

# Clustering
## 16um
se.16 <- SubsetSTData(se.hd, expression = resolution == "016um")
se.16 <- se.16 |>
    NormalizeData() |>
    ScaleData() |>
    FindVariableFeatures() |>
    RunPCA() |>
    FindNeighbors(reduction = "pca", dims = 1:10) |>
    FindClusters(resolution = 0.2)
## 8um
se.8 <- SubsetSTData(se.hd, expression = resolution == "08um")
se.8 <- se.8 |>
    NormalizeData() |>
    ScaleData() |>
    FindVariableFeatures() |>
    RunPCA() |>
    FindNeighbors(reduction = "pca", dims = 1:10) |>
    FindClusters(resolution = 0.2)

clusters_16 <- FetchData(se.16, "seurat_clusters") %>%
    rownames_to_column(var = "spot_id")
clusters_8 <- FetchData(se.8, "seurat_clusters") %>%
    rownames_to_column(var = "spot_id")
write_csv(clusters_16, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/clusters_16.csv")
write_csv(clusters_8, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/clusters_8.csv")

genes_16 <- se.16@assays$Spatial@meta.data$var.features[!is.na(se.16@assays$Spatial@meta.data$var.features)] # se.16@assays$Spatial@meta.data$var.features[!is.na(se.16@assays$Spatial@meta.data$var.features)]
genes_8 <- se.8@assays$Spatial@meta.data$var.features[!is.na(se.8@assays$Spatial@meta.data$var.features)]
write.csv(genes_16, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/genes_16.csv")
write.csv(genes_8, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/genes_8.csv")

p1 <- MapLabels(se.16,
    column_name = "seurat_clusters", override_plot_dims = TRUE,
    shape = "tile", image_use = "raw", label_by = "resolution"
)
ggsave(plot = p1, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/clusters_16.jpg", width = 10, height = 10)
p1 <- MapLabels(se.8,
    column_name = "seurat_clusters", override_plot_dims = TRUE,
    shape = "tile", image_use = "raw", label_by = "resolution"
)
ggsave(plot = p1, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/clusters_8.jpg", width = 10, height = 10)


# NMF
library(singlet)
## 16um
se.16[["Spatial3"]] <- as(
    object = subset(se.16[["Spatial"]], features = VariableFeatures(se.16)),
    Class = "Assay"
)
set.seed(42)
se.16 <- RunNMF(se.16, assay = "Spatial3", reps = 3)
## 8um
se.8[["Spatial3"]] <- as(
    object = se.8[["Spatial"]], #subset(se.8[["Spatial"]], features = VariableFeatures(se.8)),
    Class = "Assay"
)
set.seed(42)
se.8 <- RunNMF(se.8, assay = "Spatial3", reps = 3)

k.16 <- ncol(se.16@reductions$nmf@feature.loadings)
k.8 <- ncol(se.8@reductions$nmf@feature.loadings)

se.16 <- LoadImages(se.16, image_height = 2000)
p1 <- MapFeatures(se.16,
    features = paste0("NMF_", 1:k.16), label_by = "resolution",
    override_plot_dims = TRUE, image_use = "raw",
    colors = c("darkblue", "cyan", "yellow", "red", "darkred"), 
    shape = "tile"
) &
    theme(plot.title = element_blank())
ggsave(plot = p1, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/NMF_16.jpg", width = 20, height = 20)
p1 <- MapFeatures(se.16,
    features = "NMF_3", label_by = "resolution",
    image_use = "raw", override_plot_dims = TRUE,
    shape = "tile", scale_alpha = TRUE,
    colors = c("darkblue", "cyan", "yellow", "red", "darkred")
)
ggsave(plot = p1, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/NMF_16.jpg", width = 10, height = 10)

se.8 <- LoadImages(se.8, image_height = 2000)
p1 <- MapFeatures(se.8,
    features = paste0("NMF_", 1:k.8), label_by = "resolution",
    override_plot_dims = TRUE, image_use = "raw",
    colors = c("darkblue", "cyan", "yellow", "red", "darkred"), 
    shape = "tile"
) &
    theme(plot.title = element_blank())
ggsave(plot = p1, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/NMF_8.jpg", width = 20, height = 20)
p1 <- MapFeatures(se.8,
    features = "NMF_2", label_by = "resolution",
    image_use = "raw", override_plot_dims = TRUE,
    shape = "tile", scale_alpha = TRUE,
    colors = c("darkblue", "cyan", "yellow", "red", "darkred")
)
ggsave(plot = p1, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/NMF_8.jpg", width = 10, height = 10)



p1 <- MapFeatures(se.16,
    features = "NMF_3", label_by = "resolution",
    image_use = NULL, override_plot_dims = TRUE,
    shape = "tile", scale_alpha = TRUE,
    colors = c("darkblue", "cyan", "yellow", "red", "darkred")
)
ggsave(plot = p1, "/data/javier.escudero/VisiumHD/data/raw/10x_mousebrain/for_semla/duba_16.jpg", width = 10, height = 10)
