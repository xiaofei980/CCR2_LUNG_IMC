BASE = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG"

# Name of input and output directories:
INPUT_PATH = file.path(BASE, "Data/Annotation_input")
OUTPUT_PATH = file.path(BASE, "Results/Annotation_output/Distance")
library(tidyverse)
library("readr")
library("ggplot2")
library("dplyr")
library(gplots) #library for plotting data
library(tools)
# install.packages('HSAUR') ### statistical analysis package
library(HSAUR)
# install.packages('pheatmap') ### pheatmap function for k-means clustering visualisation
library('pheatmap')
# install.packages('Rtsne') ### installing packages and libraries for tSNE plots in R
library(Rtsne)
library(scales)
# install.packages('colorRamps') ### installing a package to produce plots with different colours
library(colorRamps)
# install.packages("corrplot")
library(corrplot)
# install.packages("Hmisc")
library(Hmisc)
library("RColorBrewer")
# install.packages("plotly")
library("plotly")
# install.packages("umap")
library(umap)
# install.packages("ggpubr")
library(ggpubr)
library(factoextra)



# Figure 5c
# Calculate the distance to the nearest cell of a cell type, 
# Or skip this step and load the dataset with distances
cd = celldata[,c("Location_Center_X", "Location_Center_Y", "annotation", "ROI_ID", "ExpGroup")]
names(cd)
# For cell n from cd, calculate distance to nearest cell of clusterID
cl_v = sort(unique(cd$annotation))
cd_result = data.frame()
for (ROI in unique(cd$ROI_ID)){
  cd_ROI = cd[which(cd$ROI_ID == ROI),]
  for (clusterID in cl_v) {
    cd_cl = cd_ROI[which(cd_ROI$annotation == clusterID, arr.ind = FALSE),]
    for (n in 1:nrow(cd_ROI)) {
      cd_ROI$temp[n] = sqrt(min((cd_cl[,1]- as.numeric(cd_ROI[n,1]))^2 + (cd_cl[,2]- as.numeric(cd_ROI[n,2]))^2))
    }
    column_name = paste("dist_cluster_", clusterID, sep = "")  ## next time can change to annotaiton 
    names(cd_ROI)[names(cd_ROI) == "temp"] = column_name
  }
  cd_result = rbind(cd_result, cd_ROI) 
}

celldata = cbind(celldata, cd_result[,grep("dist_cluster", names(cd_result))])

# Load celldata dataset including cluster distances
OUTPUT_PATH = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Distance/"
date<- format(Sys.Date(), "%y%m%d")
output_file = paste0(OUTPUT_PATH, date, "celldata_distance.csv")#OUTPUT_PATH
output_file
write.csv(celldata, file = output_file, row.names = F) 
celldata = read.csv(filenm)

#distance to CD8 T cells
cluster_order = c("T cell CD4","Tregs","B cells",
                  "Dendritic cells CD103", "Dendritic cells other", "Alvelor Macrophages", "Interstitial Macrophages","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour")
p = ggplot(celldata[which(!celldata$annotation %in% c("Uncertain", "T cell CD8")),], 
           aes(x = annotation, y=`dist_cluster_T cell CD8`, fill = factor(ExpGroup, levels = c("WT", "CCR2KO")))) + 
  geom_boxplot() +
  scale_color_gradientn(colours = rainbow(5)) +
  scale_x_discrete(limits = cluster_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right") +
  ylab("Distance to nearest CD8 T cell (px)") +
  scale_y_continuous(trans='log2') 
p
filename = "distance_metacluster2_to_CD8.pdf"
ggsave(plot = p, device = "pdf", width=9, height=5, dpi=300, path = OUTPUT_PATH, filename = filename)


# Heatmap of T cell neighbours
unique_clusters =  c("T cell CD4","Tregs","B cells",
                     "Dendritic cells CD103", "Dendritic cells other", "Alvelor Macrophages", "Interstitial Macrophages","Neutrophils",
                     "Endothelium","Epithelium","Fibroblasts","Tumour")

bins = list(c(0,20), c(20,40), c(40,60), c(60,80), c(80,100))
my_palette = colorRampPalette(c("yellow", "purple"))(n = 199)
for (celltype in c("T cell CD4", "T cell CD8", "Tregs")){
  dist = paste("dist_cluster_", celltype, sep = "")
  for (domain in c("Normal", "Interface", "Tumour")){
    Neighbours_matrix = data.frame()
    for (tr in c("WT", "CCR2KO")){
      for (bin in bins){
        print(bin)
        for (cl in unique_clusters){
          Neighbours_matrix[cl,paste(tr, " ", bin[1], "-", bin[2], "px", sep = "")] =  tally(celldata[which(celldata[,dist]  > bin[1] & celldata[,dist] < bin[2] &
                                                                                                              celldata$ExpGroup == tr & celldata$annotation == cl 
                                                                                                            & celldata$new_domain == domain
          ),])
        }
        
      }
    }
    Neighbours_matrix = as.matrix(Neighbours_matrix)
    pdf(file=paste(OUTPUT_PATH, celltype, "_", domain, "_Neighbours_matrix.pdf", sep = ""), width=10, height=15)
    heatmap.2(Neighbours_matrix,
              col= my_palette, scale = "col",Rowv = NA,Colv = NA, cexRow = 2.5, cexCol = 2.5,
              density.info="none",  # turns off density plot inside color legend
              trace="none", adjCol = c(1,0.5),
              colsep = c(5), sepwidth = c(0.05,0),
              key = TRUE,  keysize=0.6,
              # lmat = rbind(c(3,4),c(2,1)), #lwid = c(1.5,4), lhei = c(1.5,4,1), #This moves the key
              margins =c(25,25)
    )
    dev.off()

  }
}



