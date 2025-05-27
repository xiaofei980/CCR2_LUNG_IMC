# Post-processing
# in this code, I will include the annotation and domian information for next procedure.
library(gplots)
library(pheatmap)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(ggridges)
library(ggpubr)
library(scales)
library(clustree)
library(cowplot)
library(rlang)
library(tiff)
library(stringr)
library(purrr)



# Path for loading data and saving results 
# Project folder path:
BASE = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG"

# Name of input and output directories:
INPUT_PATH = file.path(BASE, "Data/Annotation_input")
OUTPUT_PATH = file.path(BASE, "Results/Annotation_output/")

# Load celldata excluding tSNE coordinates  /load sampled celldata with tSNE coordinates. 
celldata = read.csv(paste0(INPUT_PATH, "/250218_normalised_concatenated_scaled_data_001threshold.csv"))
#celldata_subsample = read.csv(paste0(INPUT_PATH,"/250107_celldata_subsamples.csv"))

# Add clustering and t-SNE columns to datafiles
# tSNEcsv is sampledtSNE 
tSNECSV = read.csv(file.path(INPUT_PATH, "250218_tSNE.csv"))
clustCSV = read.csv(file.path(INPUT_PATH, "250218_clustering_k30_k25_k20_k15_19markers.csv"))
celldata = cbind(celldata, clustCSV)
# Load in the image metadata
imagemetadata = read.csv(file.path(INPUT_PATH, "imagemetadata.csv"))
# Due to the computer is not powerful enough ,so I only get sampled tSNE here. celldata_subsample is using for making plots which includes tSNE coordinates.
#celldata_subsample = cbind(celldata_subsample, tSNECSV)
celldata = cbind(celldata, tSNECSV)
##celldata with part tSNE coordinates "NA" but i dont need using tSNE coordinate for heatmap
# columns_to_add <- c("cell_ID","tSNE1", "tSNE2")
# celldata <- left_join(x = celldata, y = celldata_subsample[,columns_to_add], 
#                      by = "cell_ID")

celldata <- celldata %>%
  mutate(name = paste(ROI_ID, filename, sep = "_"))
unique(celldata$name)


create.resultsdir <- function(path, name) {
  dir_path <- file.path(path, name)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
  return(dir_path)
}

date = str_sub(gsub("-", "", Sys.Date()), -6, -1)

## Due to there are still some extrame values for CD103 . so check before annotation
top20_cd103 <- celldata %>%
  arrange(desc(MI_CD103)) %>%
  slice_head(n = 20)
print(top20_cd103)



# # Manually assign cluster labels ,current K20
celldata$anncluster = stringr::str_pad(celldata$k20, 2, side="left", pad="0") 
unique(celldata$anncluster)
# Split clusters based on thresholds for marker intensities
celldata[which(celldata$anncluster == "12"), "anncluster"] <- "12B"
celldata[which(celldata$anncluster == "12B" & celldata$MI_Foxp3 >= 0.3), "anncluster"] <- "12A"
celldata[which(celldata$anncluster == "12B" & 
                 celldata$MI_CD8a >= 0.5 & celldata$MI_CD4 < 1), "anncluster"] <- "12C"

##  21  split  to neutrophils and CD11b single positive celsl


#split 21 tomono an neutropihls
celldata[which(celldata$anncluster == "21"
               & celldata$MI_LY6G > 0.3), "anncluster"] = "21A" # cehckede 0.25  similar to 0.3   0.4 is too high 
celldata[which(celldata$anncluster == "21"), "anncluster"] = "21B"


#split 17 to DC103 and uncertain -->why?For get the reason --> there are some artifacts from CD103 . 
celldata[which(celldata$anncluster == "17"
               & celldata$MI_CD103 > 200), "anncluster"] = "17A" #1.5 in previous , but incluse some cDC1
celldata[which(celldata$anncluster == "17"), "anncluster"] = "17B"



#250218 varified by xiaofei 

celldata[which(celldata$anncluster == "01"), "anncluster"] = "01_B cells" # checked 
celldata[which(celldata$anncluster == "02"), "anncluster"] = "02_Tumour" # checked
celldata[which(celldata$anncluster == "03"), "anncluster"] = "03_Tumour"#c checked
celldata[which(celldata$anncluster == "04"), "anncluster"] = "04_Tumour"# Also see loe MHCII expression ,but negative for all other markers
celldata[which(celldata$anncluster == "05"), "anncluster"] = "05_Alvelor Macrophages"#with highly CD11c, CD68+CD11c+
celldata[which(celldata$anncluster == "06"), "anncluster"] = "06_Interstitial Macrophages"# checked
celldata[which(celldata$anncluster == "07"), "anncluster"] = "07_Tumour"#checked
celldata[which(celldata$anncluster == "08"), "anncluster"] = "08_Fibroblasts"#checked
celldata[which(celldata$anncluster == "09"), "anncluster"] = "09_Dendritic cells other"# checked with febe, from umap it is dc 
celldata[which(celldata$anncluster == "10"), "anncluster"] = "10_Tumour" # umap shows more like tumour ,also not typical CK19+area 
celldata[which(celldata$anncluster == "11"), "anncluster"] = "11_Endothelium"# also express ly6c
celldata[which(celldata$anncluster == "12A"), "anncluster"] = "12A_Tregs" #
celldata[which(celldata$anncluster == "12B"), "anncluster"] = "12B_T cell CD4" #
celldata[which(celldata$anncluster == "12C"), "anncluster"] = "12C_T cell CD8" #
celldata[which(celldata$anncluster == "13"), "anncluster"] = "13_Endothelium"#Endotheliem ,also MHCII+
celldata[which(celldata$anncluster == "14"), "anncluster"] = "14_Tumour"
celldata[which(celldata$anncluster == "15"), "anncluster"] = "15_Tumour"
celldata[which(celldata$anncluster == "16"), "anncluster"] = "16_Neutrophils"# almost all express ly6g
celldata[which(celldata$anncluster == "17A"), "anncluster"] = "17A_Artifacts"
celldata[which(celldata$anncluster == "17B"), "anncluster"] = "17B_Dendritic cells CD103"
celldata[which(celldata$anncluster == "18"), "anncluster"] = "18_T cell CD8"
celldata[which(celldata$anncluster == "19"), "anncluster"] = "19_Epithelium"#mainlyCLK19 ,clear epitheliums, airway 
celldata[which(celldata$anncluster == "20"), "anncluster"] = "20_Tumour"  
celldata[which(celldata$anncluster == "21A"), "anncluster"] = "21A_Neutrophils"# will split to neutrphils and monocytes
celldata[which(celldata$anncluster == "21B"), "anncluster"] = "21B_Monocytes" # cutoff setting by images
celldata[which(celldata$anncluster == "22"), "anncluster"] = "22_Interstitial Macrophages" # mainly in one ROI _ 3256.6e.  no cD11c signal
celldata[which(celldata$anncluster == "23"), "anncluster"] = "23_Tumour" # few
celldata[which(celldata$anncluster == "24"), "anncluster"] = "24_Endothelium" # or tumour?
celldata[which(celldata$anncluster == "25"), "anncluster"] = "25_Tumour" # checked ,F480 is very low and mainly spillover
celldata[which(celldata$anncluster == "26"), "anncluster"] = "26_Endothelium"  # check cutoff 
celldata[which(celldata$anncluster == "27"), "anncluster"] = "27_Endothelium" #  also express ly6c
celldata[which(celldata$anncluster == "28"), "anncluster"] = "28_Fibroblasts" #  
celldata[which(celldata$anncluster == "29"), "anncluster"] = "29_Tumour"# no cluster in 3371 5a mainly554b-2 
celldata[which(celldata$anncluster == "30"), "anncluster"] = "30_Interstitial Macrophages"#all few checked
celldata[which(celldata$anncluster == "31"), "anncluster"] = "31_Epithelium"#no cluster in 33715a. all few
celldata[which(celldata$anncluster == "32"), "anncluster"] = "32_Uncertain" #TCRgrt?   not much in all images double positive
celldata[which(celldata$anncluster == "33"), "anncluster"] = "33_Fibroblasts"
celldata[which(celldata$anncluster == "34"), "anncluster"] = "34_Uncertain"#  mainly edge, possible real uncertain

table(celldata$anncluster)

unique(celldata$anncluster)
# Generate a column "cluster2" that makes a string out of the cluster numbers
# celldata$anncluster = annotate.cluster(celldata$anncluster)

celldata$annotation = str_split(celldata$anncluster,
                                pattern = "_", simplify = TRUE)[,2]
unique(celldata$annotation)
annlevels = c("Tumour", "Fibroblasts", "Endothelium", "Epithelium",
              "T cell CD8", "T cell CD4", "Tregs", "Dendritic cells other",
              "Dendritic cells CD103", "Alvelor Macrophages", "Interstitial Macrophages","Neutrophils", "Monocytes","B cells", "Uncertain","Artifacts")
celldata$annotation = factor(celldata$annotation, levels = annlevels)
table(celldata$anncluster)
table(celldata$annotation)
celldata_count_perROI= celldata %>%
  group_by(anncluster,ROI_ID) %>%
  dplyr::summarise(totaln=n())
View(celldata_count_perROI)


## change  few domain segmentation due to obvous error from cell profiler pipeline 

#### 32566e
normal_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/BRAC3256_6e_WT_ROI_001_1.txt_crop1/roi_brac3256_6e_wt_roi_001_1_crop1/normal_minus_interface_domain.tiff")
tumor_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/BRAC3256_6e_WT_ROI_001_1.txt_crop1/roi_brac3256_6e_wt_roi_001_1_crop1/tumour_minus_interface_domain.tiff")
interface_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/BRAC3256_6e_WT_ROI_001_1.txt_crop1/roi_brac3256_6e_wt_roi_001_1_crop1/interface_domain.tiff")
# 图像尺寸
img_height <- dim(normal_mask)[1]
img_width <- dim(normal_mask)[2]
unique(celldata$filename)
celldata_3256<- subset(celldata,filename=="ROI_BRAC3256_6E_WT_ROI1_x1")
# 可选：检查 celldata 坐标是否超出图像范围
if (any(celldata_3256$x > img_width | celldata$y > img_height)) {
  stop("celldata 中的坐标超出图像范围")
}

# Step 3: 提取每个点对应的 mask 值
celldata_3256 <- celldata_3256 %>%
  mutate(
    # TIFF 是左上角为 (1,1)，所以 y 直接用（行号），x 是列号
    normal_val = normal_mask[cbind(Location_Center_Y, Location_Center_X)],
    tumor_val = tumor_mask[cbind(Location_Center_Y, Location_Center_X)],
    interface_val = interface_mask[cbind(Location_Center_Y, Location_Center_X)],
    
    new_domain = case_when(
      interface_val == 1 ~ "Interface",
      tumor_val == 1 ~ "Tumour",
      normal_val == 1 ~ "Normal",
      TRUE ~ "n/a"
    )
  ) %>%
  select(-normal_val, -tumor_val, -interface_val)
colnames(celldata_3256)
table(celldata_3256$new_domain)

### 554b
normal_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/AWBA55_4b_CCR2KO_ROI_002_2/roi_awba55_4b_ccr2ko_roi_002_2/normal_minus_interface_domain.tiff")
tumor_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/AWBA55_4b_CCR2KO_ROI_002_2/roi_awba55_4b_ccr2ko_roi_002_2/tumour_minus_interface_domain.tiff")
interface_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/AWBA55_4b_CCR2KO_ROI_002_2/roi_awba55_4b_ccr2ko_roi_002_2/interface_domain.tiff")
# 图像尺寸
img_height <- dim(normal_mask)[1]
img_width <- dim(normal_mask)[2]
unique(celldata$filename)
celldata_554b_002<- subset(celldata,filename=="ROI_AWBA55_4B_CCR2KO_ROI2_2")
# 可选：检查 celldata 坐标是否超出图像范围
if (any(celldata_554b_002$x > img_width | celldata$y > img_height)) {
  stop("celldata 中的坐标超出图像范围")
}

# Step 3: 提取每个点对应的 mask 值
celldata_554b_002 <- celldata_554b_002 %>%
  mutate(
    # TIFF 是左上角为 (1,1)，所以 y 直接用（行号），x 是列号
    normal_val = normal_mask[cbind(Location_Center_Y, Location_Center_X)],
    tumor_val = tumor_mask[cbind(Location_Center_Y, Location_Center_X)],
    interface_val = interface_mask[cbind(Location_Center_Y, Location_Center_X)],
    
    new_domain = case_when(
      interface_val == 1 ~ "Interface",
      tumor_val == 1 ~ "Tumour",
      normal_val == 1 ~ "Normal",
      TRUE ~ "n/a"
    )
  ) %>%
  select(-normal_val, -tumor_val, -interface_val)
colnames(celldata_554b_002)
table(celldata_554b_002$new_domain)

## 3310 6a
normal_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/BRAC3310_6a_WT_ROI_001_1.txt_crop1/roi_brac3310_6a_wt_roi_001_1_crop1/normal_minus_interface_domain.tiff")
tumor_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/BRAC3310_6a_WT_ROI_001_1.txt_crop1/roi_brac3310_6a_wt_roi_001_1_crop1/tumour_minus_interface_domain.tiff")
interface_mask <- readTIFF("/mnt/Data1/groupfebe/runs/Data/XIAYU/CCR2_LUNG/results/segmentation/BRAC3310_6a_WT_ROI_001_1.txt_crop1/roi_brac3310_6a_wt_roi_001_1_crop1/interface_domain.tiff")
# 图像尺寸
img_height <- dim(normal_mask)[1]
img_width <- dim(normal_mask)[2]
unique(celldata$filename)
celldata_3310<- subset(celldata,filename=="ROI_BRAC3310_6A_WT_ROI1_x1")
# 可选：检查 celldata 坐标是否超出图像范围
if (any(celldata_3310$x > img_width | celldata$y > img_height)) {
  stop("celldata 中的坐标超出图像范围")
}

# Step 3: 提取每个点对应的 mask 值
celldata_3310 <- celldata_3310 %>%
  mutate(
    # TIFF 是左上角为 (1,1)，所以 y 直接用（行号），x 是列号
    normal_val = normal_mask[cbind(Location_Center_Y, Location_Center_X)],
    tumor_val = tumor_mask[cbind(Location_Center_Y, Location_Center_X)],
    interface_val = interface_mask[cbind(Location_Center_Y, Location_Center_X)],
    
    new_domain = case_when(
      interface_val == 1 ~ "Interface",
      tumor_val == 1 ~ "Tumour",
      normal_val == 1 ~ "Normal",
      TRUE ~ "n/a"
    )
  ) %>%
  select(-normal_val, -tumor_val, -interface_val)
colnames(celldata_3310)
table(celldata_3310$new_domain)

###### merge new domain back

# 初始化 new_domain 列
celldata$new_domain <- celldata$domain
celldata$new_domain <- as.character(celldata$domain)


# 创建一个 list 存放所有 ROI 的子集数据（注意都要包含 new_domain 列）
roi_list <- list(celldata_554b_002, celldata_3256, celldata_3310)

# 依次替换每个 ROI 的对应部分
for (roi_df in roi_list) {
  celldata <- celldata %>%
    left_join(
      roi_df %>% select(filename, Location_Center_X, Location_Center_Y, new_domain),
      by = c("filename", "Location_Center_X", "Location_Center_Y"),
      suffix = c("", ".updated")
    ) %>%
    mutate(
      new_domain = ifelse(!is.na(new_domain.updated), new_domain.updated, new_domain)
    ) %>%
    select(-new_domain.updated)
}

unique(celldata$new_domain)
table(celldata$domain)
table(celldata$new_domain)


## check from ggplot 

#original domain
p<-ggplot(celldata, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(data = celldata[which(celldata$domain=="Normal"),], size = 0.01, colour = "red") +
  geom_point(data = celldata[which(celldata$domain=="Tumour"),], size = 0.01, colour = "green") +
  geom_point(data = celldata[which(celldata$domain=="Interface"),], size = 0.01, colour = "yellow") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(),
        panel.spacing = unit(0, "mm"),
        #plot.background = element_blank()
        plot.background = element_rect(fill = "white", color = NA)
        )+
  facet_wrap(~filename,nrow=3)
print(p)
# 选择输出格式：png, tiff, 或 pdf
format <- "png"  # 修改为 "tiff" 或 "pdf" 可切换格式
filename <- paste0("/Domain_ggplot/domain_cellprofiler_plot.", format)
full_output_path <- file.path(OUTPUT_PATH, filename)

# 保存图像
ggsave(
  filename = full_output_path,
  plot = p,
  width = 10, height = 8,dpi=600
)


#new domain
ggplot(celldata, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(data = celldata[which(celldata$new_domain=="Normal"),], size = 0.01, colour = "red") +
  geom_point(data = celldata[which(celldata$new_domain=="Tumour"),], size = 0.01, colour = "green") +
  geom_point(data = celldata[which(celldata$new_domain=="Interface"),], size = 0.01, colour = "yellow") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank())+
  facet_wrap(~ROI_ID,nrow=3)

table(celldata$annotation)
celldata_clean <- subset(celldata, annotation != "Artifacts")
table(celldata_clean$annotation)
annlevels = c("Tumour", "Fibroblasts", "Endothelium", "Epithelium",
              "T cell CD8", "T cell CD4", "Tregs", "Dendritic cells other",
              "Dendritic cells CD103", "Alvelor Macrophages", "Interstitial Macrophages","Neutrophils", "Monocytes","B cells", "Uncertain")
celldata_clean$annotation = factor(celldata_clean$annotation, levels = annlevels)

date<- format(Sys.Date(), "%y%m%d")
output_file = paste0(OUTPUT_PATH,"/", date, "_celldata_clean.csv")#OUTPUT_PATH
output_file
write.csv(celldata_clean, file = output_file, row.names = F) 

table(celldata_clean$anncluster)
table(celldata_clean$annotation)



