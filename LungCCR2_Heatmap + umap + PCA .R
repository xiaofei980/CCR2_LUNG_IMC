##Heatmap
library(gplots)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(umap)

BASE = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG"

# Name of input and output directories:
INPUT_PATH = file.path(BASE, "Data/Annotation_input")
OUTPUT_PATH = file.path(BASE, "Results/Annotation_output/Heatmap/")
celldata = read.csv(paste0(INPUT_PATH, "/250526_celldata_clean.csv"))

# In this step, the order is also important 
lin_markers = c("CD45", "CD3e", "CD4", "Foxp3", "CD8a","B220",  "CD11c", "MHCII","CD103", 
                "F480","CD68", "CD11b", "LY6G","Ly6C","EpCAM","CK19","aSMA","PECAM","CD44")
#checking the name


excl_cols = c("MI_TumourminusInterface", "MI_NormalminusInterface","MI_StructuralImage",
              "MI_InterfaceImage", "MI_Xenon131", "MI_Xenon132", "MI_DNA1",
              "MI_DNA2", "MI_Argon80")
mismatch_names <- setdiff(excl_cols, colnames(celldata))
print(mismatch_names) 
all_markers = grep("MI_", names(celldata), value = T)
all_markers = all_markers[!all_markers %in% excl_cols] %>%
  str_remove(pattern = "MI_")

lin_markers_check<- lin_markers[!lin_markers %in% all_markers]
lin_markers_check

##defien heatmap based on Complex heatmap 
plot.heatmap_CH <- function(celldata, markerlist, cluster_col,
                            title = "Marker Expression Heatmap",
                            Colv = TRUE, Rowv = TRUE,
                            transpose = FALSE,
                            scale_by = "marker",  # NEW: "marker", "cluster", or FALSE
                            save_path = NULL,
                            file_type = c("png", "jpeg", "pdf", "tiff", "svg"),
                            width = 10, height = 7, dpi = 300,
                            also_plot = TRUE) {
  
  file_type <- match.arg(file_type)
  
  # 1. Remove rows with missing cluster info
  celldata <- celldata[!is.na(celldata[[cluster_col]]), ]
  
  # 2. Calculate mean expression of each marker per cluster
  cluster_means <- celldata %>%
    group_by(!!sym(cluster_col)) %>%
    summarise(across(all_of(paste0("MI_", markerlist)), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup()
  
  # 3. Convert to matrix
  mat <- as.data.frame(cluster_means)
  rownames(mat) <- mat[[cluster_col]]
  mat[[cluster_col]] <- NULL
  mat <- as.matrix(mat)
  
  # 4. Clean column names
  colnames(mat) <- gsub("^MI_", "", colnames(mat))
  
  # 5. Force marker column order
  mat <- mat[, markerlist, drop = FALSE]
  
  # 6. Scale logic
  if (is.character(scale_by)) {
    if (scale_by == "marker") {
      mat <- scale(mat)
    } else if (scale_by == "cluster") {
      mat <- t(scale(t(mat)))
    } else {
      stop("scale_by must be 'marker', 'cluster', or FALSE")
    }
  }
  
  # 7. Transpose if needed
  if (transpose) {
    mat <- t(mat)
  }
  
  # 8. Heatmap settings
  cluster_rows <- Rowv
  cluster_columns <- Colv
  
  row_angle <- if (transpose) 45 else 0
  row_fontsize <- if (transpose) 8 else 10
  row_side <- if (transpose) "right" else "left"
  
  ht <- Heatmap(mat,
                name = "Expr",
                col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                cluster_rows = cluster_rows,
                cluster_columns = cluster_columns,
                column_names_rot = if (transpose) 45 else 90,
                column_title = title,
                row_title = NULL,
                heatmap_legend_param = list(title = "Z-score"),
                column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                column_names_gp = gpar(fontsize = 10),
                row_names_gp = gpar(fontsize = row_fontsize, rot = row_angle),
                row_names_side = row_side)
  
  # 9. Save to file
  if (!is.null(save_path)) {
    switch(file_type,
           png  = png(paste0(save_path, ".png"), width = width, height = height, units = "in", res = dpi),
           jpeg = jpeg(paste0(save_path, ".jpeg"), width = width, height = height, units = "in", res = dpi),
           pdf  = pdf(paste0(save_path, ".pdf"), width = width, height = height),
           tiff = tiff(paste0(save_path, ".tiff"), width = width, height = height, units = "in", res = dpi),
           svg  = svg(paste0(save_path, ".svg"), width = width, height = height)
    )
    draw(ht, padding = unit(c(2, 2, 2, 6), "cm"))
    dev.off()
  }
  
  # 10. Optional viewer display
  if (also_plot) {
    draw(ht, padding = unit(c(2, 2, 2, 6), "cm"))
  }
}
plot.heatmap_CH(
  celldata = celldata,
  markerlist = lin_markers,
  cluster_col = "anncluster",
  title = "Lineage Marker Heatmap",
  scale_by = "cluster",
  Colv = FALSE,
  Rowv = FALSE,
  transpose = TRUE,
  save_path = file.path(OUTPUT_PATH, "lineage_anncluster_heatmap_transport"),
  file_type = "svg"
)
  
plot.heatmap_CH(
  celldata = celldata,
  markerlist = lin_markers,
  cluster_col = "annotation",
  title = "Lineage Marker Heatmap",
  scale_by = "cluster",
  Colv = FALSE,
  Rowv = FALSE,
  transpose = TRUE,
  save_path = file.path(OUTPUT_PATH, "lineage_annotation_heatmap_transport"),
  file_type = "svg"
)


#_UMAP_______________________________-

library(gplots)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(gridExtra)
library(ggridges)
library(ggpubr)
library(scales)
library(cowplot)
library(umap)
library(magick)

OUTPUT_PATH = file.path(BASE, "Results/Annotation_output/UMAP/")

generate.UMAP <- function(celldata, clusters, markers, nsample, filename) {
  # Filter and sample data
  data <- celldata[celldata$anncluster %in% clusters, ]
  data <- data[sample(x = nrow(data), size = nsample), ]
  
  # Run UMAP
  umapResults <- umap(data[, markers], verbose = TRUE, n_neighbors = 10)
  data$umap1 <- umapResults$layout[, 1]
  data$umap2 <- umapResults$layout[, 2]
  
  # Save with date-prefixed filename
  dated_filename <- file.path(OUTPUT_PATH, paste0(Sys.Date(), "_", filename))
  write.csv(x = data, file = dated_filename, row.names = TRUE)
  
  return(data)
}


# in case there is no get clusters function defined.
get.clusters = function(patterns) {
  clusters = sapply(X = patterns, FUN = grep, x = unique(celldata$anncluster), value = T)
  clusters = unlist(clusters, use.names = F)
  return(clusters)
} 


### call the define umap function 

TCellClusters = get.clusters(patterns = c("T cell","Tregs"))
TCellMarkers_lin = c("MI_CD3e", "MI_CD45", "MI_CD4", "MI_CD8a", "MI_Foxp3", "MI_TCRgdt")
TCellMarkers = c("MI_CD3e", "MI_CD45", "MI_CD4", "MI_CD8a", "MI_Foxp3",
                 "MI_LAG3", "MI_PD1", "MI_TCRgdt", "MI_TIM3")

print(TCellClusters)
set.seed(222)
TCD = generate.UMAP(celldata = celldata, clusters = TCellClusters,
                    markers = TCellMarkers_lin, nsample = 5000, filename = "TCD.csv")


set.seed(111)  
# Macrophages & Dendritic cells
MacDendNeuMonoCellClusters = get.clusters(patterns = c("Macrophage", "Dendritic","DC","Neutrophils","Mono"))
MacDendNeuMonoCell_lineage = c("MI_MHCII", "MI_F480", "MI_CD103", "MI_CD44",
                               "MI_CD68", "MI_CD206",  "MI_CD11c", "MI_CD11b","MI_LY6G") 
MacDendNeuMonoCellColumns = c("MI_MHCII", "MI_F480", "MI_CD103", "MI_PDL1", "MI_CD86","MI_Ly6C",
                              "MI_CD68", "MI_CD206",  "MI_CD11c", "MI_CD11b","MI_LY6G",  "MI_CD44","MI_Sirpa") 
MacDendNeuMonoCD = generate.UMAP(celldata = celldata, clusters = MacDendNeuMonoCellClusters,
                                 markers = MacDendNeuMonoCell_lineage, nsample = 20000,
                                 filename = "MacDendNeuMonoCD.csv")

# Macrophages & Dendritic cells
MacDendCellClusters = get.clusters(patterns = c("Macrophage", "Dendritic"))
print(MacDendCellClusters)
MacDendCellColumns_lin= c("MI_MHCII", "MI_F480", "MI_CD103", 
                          "MI_CD68", "MI_CD11c", "MI_CD11b")
MacDendCellColumns = c("MI_MHCII", "MI_F480", "MI_CD103", "MI_PDL1", "MI_CD86",
                       "MI_CD68", "MI_CD206",  "MI_CD11c", "MI_CD11b","MI_LY6G","MI_CD44","MI_Sirpa") #, "MI_LY6G",  "MI_CD44")
set.seed(222)
MacDendCD = generate.UMAP(celldata = celldata, clusters = MacDendCellClusters,
                          markers = MacDendCellColumns_lin, nsample = 20000,
                          filename = "MacDendCD.csv")


MacCellClusters = get.clusters(patterns = c("Macrophage"))
print(MacCellClusters)
MacCellColumns_lin= c("MI_MHCII", "MI_F480", "MI_CD206",
                          "MI_CD68", "MI_CD11c", "MI_CD11b")
MacCellColumns = c("MI_MHCII", "MI_F480",  "MI_PDL1", "MI_CD86",
                       "MI_CD68", "MI_CD206",  "MI_CD11c", "MI_CD11b","MI_CD44","MI_Sirpa") #, "MI_LY6G",  "MI_CD44")
set.seed(222)
MacCD = generate.UMAP(celldata = celldata, clusters = MacCellClusters,
                          markers = MacCellColumns_lin, nsample = 20000,
                          filename = "MacCD.csv")


set.seed(222)
AllClusters = get.clusters(patterns = c("Neutrophils", "Monocytes", "Fibroblasts", "Macrophages", 
                                        "Dendritic", "T cell", "Tregs", "Tumour", 
                                        "Epithelium", "Endothelium", "Uncertain"))
lineage_marker=c("MI_CD45", "MI_CD4", "MI_CD3e", "MI_CD8a","MI_B220","MI_Foxp3",
                 "MI_CD68",  "MI_EpCAM", "MI_F480",  "MI_LY6G","MI_CD103", "MI_CD11b", "MI_CD11c", "MI_MHCII", "MI_Ly6C",
                 "MI_PECAM", "MI_aSMA",  "MI_CD44", "MI_CK19")
CellColumns<-lineage_marker
excl_cols = c("MI_TumourminusInterface", "MI_NormalminusInterface","MI_StructuralImage",
              "MI_InterfaceImage", "MI_Xenon131", "MI_Xenon132", "MI_DNA1",
              "MI_DNA2", "MI_Argon80")
mismatch_names <- setdiff(excl_cols, colnames(celldata))
print(mismatch_names) 
all_markers = grep("MI_", names(celldata), value = T)
all_markers = all_markers[!all_markers %in% excl_cols] 
CellCD =generate.UMAP(celldata = celldata, clusters = AllClusters,
                      markers = lineage_marker, nsample = 20000,
                      filename = "CellCD.csv")






# Define the function to plot UMAP
plot_umap <- function(data, color_var, plot_title, legend_title, file_name, subfolder) {
  # Create the ggplot
  p <- ggplot(data, aes(x=umap1, y=umap2, colour = as.factor(data[[color_var]]))) + 
    geom_point(size = 0.3) + 
    ggtitle(plot_title) + 
    guides(colour = guide_legend(override.aes = list(size = 5), title = legend_title)) +
    theme_classic() + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.background = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          panel.grid = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_manual(values = cluster_cols)
  # Save the plot
  ggsave(filename = file_name,
         path = file.path(OUTPUT_PATH,subfolder),
         bg = "white",
         width = 2100, height = 1700, units = "px")
  return(p)
}

cluster_cols = c( "#FFF4A4FF",  "#0072B2", "#D55E00", "#CC79A7", "#ABDDA4FF","#6A5ACD", "#618F75FF",  "#DC143C",  "#40E0D0","#336666",
                  "#4E79A7FF", "#A0CBE8FF","#938ABBFF",  "#B07AA1FF",
                  "#FF9D9AFF", "#CC6666",  "#704850FF",  "#8C8C8C", "#FFD700",   "#00FF7F",   "#8B4513",  "#FF6347",    "#7FFF00",  
                  "#8B0000","#F0E442", 'pink','purple','blue',"#8FBC8B",  "#4682B4", "#FF4500", "#2E8B57",
                  "#483D8B", "#FF1493", "#FF8C00", "#98FB98", "#556B2F",
                  "#CD5C5C", "#00008B", "#FFDAB9", "#BA55D3", "#87CEEB",
                  "#B0E0E6", "#5F9EA0")

# T cell 
TCD_filter<-TCD
TCD_filter= TCD_filter[TCD_filter$umap1<10,]
TCD_filter= TCD_filter[TCD_filter$umap2<15,]
plot_umap(TCD_filter,color_var = "anncluster",plot_title="U-MAP T cells",legend_title = "Clusters",file_name="UMAP_T_cell_cluster.jpeg",subfolder = "UMAP/T_Cell_markers/Montage/")
plot_umap(TCD_filter,color_var = "annotation",plot_title="U-MAP T cells",legend_title = "Annotation",file_name="UMAP_T_cell_annotation.jpeg",subfolder = "UMAP/T_Cell_markers/Montage/")
plot_umap(TCD_filter,color_var = "ExpGroup",plot_title="U-MAP T cells",legend_title = "Group",file_name="UMAP_T_cell_group.jpeg",subfolder = "UMAP/T_Cell_markers/Montage/")

# Myloid cell
MacCD_filter<-MacCD
MacCD_filter= MacCD_filter[MacCD_filter$umap1< 15,]
MacCD_filter= MacCD_filter[MacCD_filter$umap1 > -15,]
MacCD_filter= MacCD_filter[MacCD_filter$umap2< 10,]
#ggplot(MacCD_filter,aes(x=umap1,y=umap2,colour="ExpGroup"))
plot_umap(MacCD_filter,color_var = "anncluster",plot_title="U-MAP Mac cells",legend_title = "Clusters",file_name="UMAP_MacCD_cluster.jpeg",subfolder = "UMAP/MacCD_markers/Montage/")
plot_umap(MacCD_filter,color_var = "annotation",plot_title="U-MAP Mac cells",legend_title = "Annotation",file_name="UMAP_MacCD_annotation.jpeg",subfolder = "UMAP/MacCD_markers/Montage/")
plot_umap(MacCD_filter,color_var = "ExpGroup",plot_title="U-MAP Mac cells",legend_title = "Group",file_name="UMAP_MacCD_group.jpeg",subfolder = "UMAP/MacCD_markers/Montage/")

# Myloid cell
MacDendCD_filter<-MacDendCD
MacDendCD_filter= MacDendCD_filter[MacDendCD_filter$umap1>-15,]
MacDendCD_filter= MacDendCD_filter[MacDendCD_filter$umap2>-15,]
#ggplot(MacDendCD_filter,aes(x=umap1,y=umap2,colour="ExpGroup"))
plot_umap(MacDendCD_filter,color_var = "anncluster",plot_title="U-MAP MacDend cells",legend_title = "Clusters",file_name="UMAP_MacDendCD_cluster.jpeg",subfolder = "UMAP/MacDendCD_markers/Montage/")
plot_umap(MacDendCD_filter,color_var = "annotation",plot_title="U-MAP MacDend cells",legend_title = "Annotation",file_name="UMAP_MacDendCD_annotation.jpeg",subfolder = "UMAP/MacDendCD_markers/Montage/")
plot_umap(MacDendCD_filter,color_var = "ExpGroup",plot_title="U-MAP MacDend cells",legend_title = "Group",file_name="UMAP_MacDendCD_group.jpeg",subfolder = "UMAP/MacDendCD_markers/Montage/")
#MAcDENd NEu +Mono
MacDendNeuMonoCD_filter<-MacDendNeuMonoCD
MacDendNeuMonoCD_filter= MacDendNeuMonoCD_filter[MacDendNeuMonoCD_filter$umap1>-15,]
MacDendNeuMonoCD_filter= MacDendNeuMonoCD_filter[MacDendNeuMonoCD_filter$umap2> -15,]
#ggplot(MacDendCD_filter,aes(x=umap1,y=umap2,colour="ExpGroup"))
plot_umap(MacDendNeuMonoCD_filter,color_var = "anncluster",plot_title="U-MAP MacDendNeu cells",legend_title = "Clusters",file_name="UMAP_MacDendNeuMonoCD_cluster.jpeg",subfolder = "UMAP/MacDendNeuMonoCD_markers/Montage/")
plot_umap(MacDendNeuMonoCD_filter,color_var = "annotation",plot_title="U-MAP MacDendNeu cells",legend_title = "Annotation",file_name="UMAP_MacDendNeuMonoCD_annotation.jpeg",subfolder = "UMAP/MacDendNeuMonoCD_markers/Montage/")
plot_umap(MacDendNeuMonoCD_filter,color_var = "ExpGroup",plot_title="U-MAP MacDendNeu cells",legend_title = "Group",file_name="UMAP_MacDendNeuMonoCD_group.jpeg",subfolder = "UMAP/MacDendNeuMonoCD_markers/Montage/")
plot_umap(MacDendNeuMonoCD_filter,color_var = "ROI_ID",plot_title="U-MAP MacDendNeu cells",legend_title = "Group",file_name="UMAP_MacDendNeuMonoCD_ROI.jpeg",subfolder = "UMAP/MacDendNeuMonoCD_markers/Montage/")

#  all
CellCD_filter<-CellCD
CellCD_filter= CellCD_filter[CellCD_filter$umap1<15,]
CellCD_filter= CellCD_filter[CellCD_filter$umap2<10,]

plot_umap(CellCD_filter,color_var = "anncluster",plot_title="U-MAP CellCD cells",legend_title = "Clusters",file_name="UMAP_CellCD_cluster.jpeg",subfolder = "UMAP/CellCDmarkers/Montage/")
plot_umap(CellCD_filter,color_var = "annotation",plot_title="U-MAP CellCD cells",legend_title = "Annotation",file_name="UMAP_CellCD_annotation.jpeg",subfolder = "UMAP/CellCDmarkers/Montage/")
plot_umap(CellCD_filter,color_var = "ExpGroup",plot_title="U-MAP CellCD cells",legend_title = "Group",file_name="UMAP_CellCD_group.jpeg",subfolder = "UMAP/CellCDmarkers/Montage/")



# plot marker expression on Umap

for (i in TCellMarkers){
  marker = str_remove(i, "MI_")
  p = ggplot(TCD_filter,
             aes(x=umap1, y=umap2, colour = get(i))) + 
    geom_point(size = 0.3) + ggtitle(paste0("T cell ", marker, " expression")) + 
    theme_classic() + labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.background = element_blank(),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          panel.grid = element_blank(),
          #axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,2))
  ggsave(filename = paste0("UMAP_T_cell_", marker, "_expression.jpeg"),
         path = file.path(OUTPUT_PATH, "UMAP/T_Cell_markers"), bg = "white",
         width = 2100, height = 1700, units = "px")
  print(p)
}
#Mac
for (i in MacCellColumns){
  marker = str_remove(i, "MI_")
  p = ggplot(MacCD_filter,
             aes(x=umap1, y=umap2, colour = get(i))) + 
    geom_point(size = 0.3) + ggtitle(paste0("MacDend cell ", marker, " expression")) + 
    theme_classic() + labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.background = element_blank(),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          panel.grid = element_blank(),
          #axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,2))
  ggsave(filename = paste0("UMAP_Mac_cell_", marker, "_expression.jpeg"),
         path = file.path(OUTPUT_PATH, "UMAP/MacCD_markers"), bg = "white",
         width = 2100, height = 1700, units = "px")
  print(p)
}
##MacDend markers
for (i in MacDendCellColumns){
  marker = str_remove(i, "MI_")
  p = ggplot(MacDendCD_filter,
             aes(x=umap1, y=umap2, colour = get(i))) + 
    geom_point(size = 0.3) + ggtitle(paste0("MacDend cell ", marker, " expression")) + 
    theme_classic() + labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.background = element_blank(),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          panel.grid = element_blank(),
          #axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,2))
  ggsave(filename = paste0("UMAP_MacDend_cell_", marker, "_expression.jpeg"),
         path = file.path(OUTPUT_PATH, "UMAP/MacDendCD_markers"), bg = "white",
         width = 2100, height = 1700, units = "px")
  print(p)
}

# MascDendMoononeu
for (i in MacDendNeuMonoCellColumns){
  marker = str_remove(i, "MI_")
  p = ggplot(MacDendNeuMonoCD_filter,
             aes(x=umap1, y=umap2, colour = get(i))) + 
    geom_point(size = 0.3) + ggtitle(paste0("MacDendNeuMono cell ", marker, " expression")) + 
    theme_classic() + labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.background = element_blank(),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          panel.grid = element_blank(),
          #axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,2))
  ggsave(filename = paste0("UMAP_MacDendNeuMono_cell_", marker, "_expression.jpeg"),
         path = file.path(OUTPUT_PATH, "UMAP/MacDendNeuMonoCD_markers"), bg = "white",
         width = 2100, height = 1700, units = "px")
  print(p)
}

#all
for (i in all_markers){
  marker = str_remove(i, "MI_")
  p = ggplot(CellCD_filter,
             aes(x=umap1, y=umap2, colour = get(i))) + 
    geom_point(size = 0.3) + ggtitle(paste0("All cell ", marker, " expression")) + 
    theme_classic() + labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.background = element_blank(),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          panel.grid = element_blank(),
          #axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,2))
  ggsave(filename = paste0("UMAP_all_cell_", marker, "_expression.jpeg"),
         path = file.path(OUTPUT_PATH, "UMAP/CellCDmarkers"), bg = "white",
         width = 2100, height = 1700, units = "px")
  print(p)
}


#-----------------------------------------------------------------------------------------------------

#Splicing the images

library(magick)

gc()  # 

# Configuration list: contains the folder, number of rows, and number of columns for each batch.
configs <- list(
  # list(
  #   input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendCD_markers/Montage/",
  #   output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendCD_markers/Montage/",
  #   n_rows = 1,
  #   n_cols = 3
  # ),
  # list(
  #   input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendCD_markers/",
  #   output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendCD_markers/",
  #   n_rows = 3,
  #   n_cols = 4
  # ),
  list(
    input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendNeuMonoCD_markers/Montage/",
    output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendNeuMonoCD_markers/Montage/",
    n_rows = 1,
    n_cols = 3
  ),
  list(
    input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendNeuMonoCD_markers/",
    output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/MacDendNeuMonoCD_markers/",
    n_rows = 3,
    n_cols = 4
  )

  
  # list(
  #   input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/CellCDmarkers/Montage/",
  #   output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/CellCDmarkers/Montage/",
  #   n_rows = 1,
  #   n_cols = 3
  # ),
  # list(
  #   input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/CellCDmarkers/1/",
  #   output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/CellCDmarkers/1/",
  #   n_rows = 3,
  #   n_cols = 4
  # )
  # list(
  #   input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/CellCDmarkers/",
  #   output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/CellCDmarkers/",
  #   n_rows = 3,
  #   n_cols = 4
  # )
  # ,
  # list(
  #   input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/T_Cell_markers/Montage/",
  #   output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/T_Cell_markers/Montage/",
  #   n_rows = 1,
  #   n_cols = 3
  # ),
  # list(
  #   input_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/T_Cell_markers/",
  #   output_folder = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/UMAP/T_Cell_markers/",
  #   n_rows = 3,
  #   n_cols = 4
  # )
  
)

today_str <- format(Sys.Date(), "%Y-%m-%d")  # get date

for (config in configs) {
  input_folder <- config$input_folder
  output_folder <- config$output_folder
  n_rows <- config$n_rows
  n_cols <- config$n_cols
  
  if (!dir.exists(input_folder)) {
    cat("Input folder does not exist:", input_folder, "\n")
    next
  }
  
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    cat("Created output folder:", output_folder, "\n")
  }
  
  jpeg_files <- list.files(input_folder, pattern = "\\.jpeg$", full.names = TRUE)
  cat("Found", length(jpeg_files), "JPEG files in", input_folder, "\n")
  
  if (length(jpeg_files) == 0) {
    next
  }
  
  images_per_combined <- n_rows * n_cols
  num_batches <- ceiling(length(jpeg_files) / images_per_combined)
  
  for (i in 1:num_batches) {
    start_idx <- (i - 1) * images_per_combined + 1
    end_idx <- min(i * images_per_combined, length(jpeg_files))
    
    batch_files <- jpeg_files[start_idx:end_idx]
    images <- lapply(batch_files, image_read)
    
    combined_image <- image_montage(
      image_join(images),
      tile = paste0(n_cols, "x", n_rows),
      geometry = "x500",
      bg = "white"
    )
    
    # add date
    output_path <- file.path(output_folder, paste0("combined_image_", today_str, "_", i, ".jpg"))
    
    image_write(combined_image, path = output_path, format = "jpg")
    cat("Image saved at:", output_path, "\n")
  }
}

##PCA------------------------------------------
library(ggplot2)
library(factoextra)  

OUTPUT_PATH = file.path(BASE, "Results/Annotation_output/PCA/")
 ##  PCA--> EPCAM plays most influence ,so here I did not include epcam 
# PCA Principal component analysis, per mouse
selectcolumns = c("ROI_ID", "ExpGroup", "MouseID",  "MI_CD45", "MI_CD3e", "MI_CD4", "MI_Foxp3", "MI_CD8a", "MI_B220",
                  "MI_CD11c", "MI_MHCII", "MI_CD103", "MI_F480", "MI_CD68", "MI_CD11b",
                  "MI_LY6G", "MI_Ly6C", "MI_EpCAM", "MI_CK19", "MI_aSMA", "MI_PECAM", "MI_CD44")
#without EPCAM
selectcolumns = c("ROI_ID", "ExpGroup", "MouseID",  "MI_CD45", "MI_CD3e", "MI_CD4", "MI_Foxp3", "MI_CD8a", "MI_B220",
                  "MI_CD11c", "MI_MHCII", "MI_CD103", "MI_F480", "MI_CD68", "MI_CD11b",
                  "MI_LY6G", "MI_Ly6C",  "MI_CK19", "MI_aSMA", "MI_PECAM", "MI_CD44")
missing_markers <- selectcolumns[!selectcolumns %in% colnames(celldata)]
if (length(missing_markers) == 0) {
  cat("All markers are present in the data frame.\n")
} else {
  cat("The following markers are missing from the data frame:\n")
  print(missing_markers)
}

cd = celldata[,selectcolumns]
PCA_input = data.frame()
MI_Ch = grep("MI_", selectcolumns, value = TRUE)
for (n in MI_Ch){
  print(n)
  for (rn in unique(cd$MouseID)){
    md = cd[which(cd$MouseID == rn),]
    PCA_input[rn,substr(n, start = 4, stop = nchar(n))] = mean(md[[n]])
    PCA_input[rn,"MouseID"] = unique(md$MouseID)
    PCA_input[rn,"ExpGroup"] = unique(md$ExpGroup)
  }
}
MI_Ch = sapply(X = MI_Ch, FUN = substr, start = 4, stop = 20)
PCA_results = prcomp(PCA_input[,MI_Ch])

p = fviz_pca_var(PCA_results,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)
print(p)
filename = "PCA_variables_lin markers.pdf"
ggsave(plot = p, device = "pdf", width=5, height=4, dpi=300, path = OUTPUT_PATH, filename = filename)


