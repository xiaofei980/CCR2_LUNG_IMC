---
title: "250407HNSCC_communityFigures"
author: "Xiaofei YU"
date: "2025-04-08"
output: html_document
---

library(ggplot2) 
library(dplyr)
library(Rtsne) 
library(clustree) 
library(ggdendro) 
library(tiff) 
library(reshape2) 
library(gplots)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(gridExtra)
library(ggridges)
library(ggpubr)
library(scales)
library(clustree)
library(cowplot)
library(umap)
library(lme4)
library(Rphenograph)
library(ggplot2)
library(dplyr) 
library(RColorBrewer) 
library(tibble)

library(devtools)
#devtools::install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap) 
library(gplots) 
library(ggpubr) 
library(textshape)
library(tibble)
library(tidyr) 


#global setting
BASE = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG"

# Name of input and output directories:
INPUT_PATH = file.path(BASE, "Data/Community_input/")
OUTPUT_PATH = file.path(BASE, "Results/Community_output")

# Name of input and output directories:
data1 = read.csv(paste0(INPUT_PATH, "/250526celldata_distance.csv"))
Rphenograph20_250523 = read.csv(paste0(INPUT_PATH, "20250526Rphenograph_CCR2LUNG_output_45clusters_k300_15ct_fractions.csv"),check.names=F)
fre_agg20 = read.csv(paste0(INPUT_PATH, "20250526_average_neighbours_45clusters_agglom_average.csv"),check.names=F)
unique(names(Rphenograph20_250523))
unique(names(fre_agg20))
unique(data1$annotation)
unique(names(data1))

#250523  new code  
# 移除没有列名的列
Rphenograph_clean <- Rphenograph20_250523[, names(Rphenograph20_250523) != ""]
fre_agg_clean <- fre_agg20

# 检查重复列（除了"cluster"）
common_cols <- intersect(names(Rphenograph_clean), names(fre_agg_clean))
common_cols <- setdiff(common_cols, "cluster")  # 保留cluster用于merge

# 从fre_agg_clean中移除重复列
fre_agg_clean <- fre_agg_clean[, !names(fre_agg_clean) %in% common_cols]

# 合并两个数据框（按cluster）
merged_data <- merge(Rphenograph_clean, fre_agg_clean, by = "cluster", all = TRUE)

# 输出合并后的结果
names(merged_data)

data1$source_ID<-data1$cell_ID
# 正确匹配列名，如果列名是 "source_ID" 而不是 "sourceID"
columns_to_add <- c("source_ID",  "new_domain", "MouseID", "ExpGroup")
celldata_subset <- data1[, columns_to_add]

# 合并时也使用 "source_ID"
merged_data <- merge(merged_data, celldata_subset, by = "source_ID", all.x = TRUE)
names(merged_data)


date<- format(Sys.Date(), "%y%m%d")
output_file = paste0(INPUT_PATH,"/", date, "_celldata_withcommunity.csv")
output_file
write.csv(merged_data, file = output_file, row.names = F) 






# barplot ----------------------------------------------------------------------------------------------


colours =  c("B cells" = "#945931",  "Dendritic cells other" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
             "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF",  
             "Alvelor Macrophages" = "#4E79A7FF", "Interstitial Macrophages" = "#5F9EA0",  "Neutrophils" = "#A0CBE8FF",  "T cell CD4" = "#B07AA1FF","Monocytes"= "#40E0D0",  "T cell CD8" = "#FF9D9AFF", "Tregs" = "#CC6666", "Tumour" = "#704850FF", "Uncertain" = "#8C8C8C")

#"Macrophages F480" = "#4E79A7FF", "Macrophages CD68" = "#5F9EA0", "Macrophages CD11b" = "#0072B2",



#ordered by T cell CD8 count, show percentage
p = ggplot(merged_data, aes(
  x = reorder(factor(agglomerate_20), source_cluster == "T cell CD8", sum),
  fill = as.factor(source_cluster)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  scale_x_discrete(limits = rev) +
  xlab("Community") +
  labs(fill = "") +
  geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution") +
  scale_fill_manual(values = colours) 
p
##ordered by T cell cD8 count ,show count 
p = ggplot(merged_data, aes(
  x = reorder(factor(agglomerate_20), source_cluster == "T cell CD8", sum),
  fill = as.factor(source_cluster)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  scale_x_discrete(limits = rev) +
  xlab("Community") +
  labs(fill = "") +
  geom_bar(stat = "count") +  # Change position to count
  ylab("Count") +  # Change y-axis label to "Count"
  scale_fill_manual(values = colours) 
p

# ordered by T cell percentage
library(dplyr)

# Calculate the percentage of T cell CD8 in each community first
# First calculate the percentage of T cell CD8 in each community
# Calculate T cell CD8 percentage for each community
community_tcd8_percentage <- merged_data %>%
  group_by(agglomerate_20) %>%
  dplyr::summarise(
    total_cells = dplyr::n(),
    tcd8_count = sum(source_cluster == "T cell CD8"),
    tcd8_percentage = tcd8_count / total_cells * 100,
    .groups = 'drop'
  ) %>%
  arrange(desc(tcd8_percentage))  # Sort by highest % first
# Add the percentage column back to the original data
merged_data <- merged_data %>%
  left_join(community_tcd8_percentage, by = "agglomerate_20")
# Then create the plot ordered by T cell CD8 percentage
p = ggplot(merged_data, aes(
  x = reorder(factor(agglomerate_20), tcd8_percentage),
  fill = as.factor(source_cluster)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  scale_x_discrete(limits = rev) +
  xlab("Community") +
  labs(fill = "") +
  geom_bar(position = "fill") + 
  scale_y_continuous(labels = scales::percent) + 
  ylab("Percentage distribution") +
  scale_fill_manual(values = colours)
p


## reorder based on th enumber of WT cells
community_counts <- fre_community20 %>%
  filter(ExpGroup == "WT") %>%
  group_by(community) %>%
  dplyr::summarise(WT_count = dplyr::n(), .groups = 'drop') %>%
  arrange(desc(WT_count))

# 将 `community` 列的顺序设置为按 `WT_count` 排序的顺序
fre_community20$community <- factor(fre_community20$community, levels = community_counts$community)
p = ggplot(fre_community20, aes(
  x = community,
  fill = as.factor(source_cluster)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16)
  ) +
  xlab("Community") +
   labs(fill = "") +
  geom_bar(stat = "count") +  # Change position to count
  ylab("Count") +  # Change y-axis label to "Count"
  scale_fill_manual(values = colours) 
p


#ordered by Group----------------------------------
fre_community20<-merged_data

#ordered by Group
# Stacked bar function 
treat_col = c("WT" = "#F8766D", "CCR2KO" = "#00BFC4")
domain_col = c( "interface" = "#00BA38", "Normal" = "#619CFF","Tumour"= "#F95")
fre_community20 <- fre_community20 %>%
  mutate(across(c("ExpGroup","annotation","new_domain"), as.factor))
str(fre_community20)
stacked_bar = function(data,
                       X,
                       fill,
                       select,
                       path,
                       ROI_name,
                       position = "fill",
                       colour,
                       width = 10,
                       height = 10,
                       split = 'single') {
  
  # 创建基础 ggplot 图形
  p = ggplot(data, aes(
    x = reorder(get(X), get(fill) == select, FUN = mean),#按每个X中指定类别（select）的平均值来排序
    fill = as.factor(get(fill))
  )) +
    theme_classic() +
    coord_flip() +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text = element_text(size = 16)
    ) +
    theme(legend.text = element_text(size = 16)) +
    xlab("Community") +
    labs(fill = "")
  
  # 根据 position 参数选择显示比例或计数
  if (position == "fill") {
    p = p + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution")
  } else if (position == "count") {
    p = p + geom_bar(stat = "count") + ylab("Cell count")
  }
  
  # 根据 split 参数决定是否拆分面板
  if (split == "multiple"){
    p = p + facet_wrap(as.formula(paste("~", fill)))  # 使用动态的列名
    p = p + theme(strip.text = element_text(size = 16))
  }
  
  # 显式设置颜色
  p = p + scale_fill_manual(values = colour)
  
  # 保存图形
  ggsave(
    plot = p,
    device = "png",
    width = width,
    height = height,
    dpi = 300,
    path = path,
    filename = paste(ROI_name, ".png", sep = "")
  )
  
  return(p)
}

# 调用函数时，你可以指定 split 和 fill 作为动态输入
stacked_bar(fre_community20, 
           X= "agglomerate_20", 
           fill= "ExpGroup", 
           select =  "WT", 
            OUTPUT_PATH, 
            paste0(date, "_Rphenograph20_250526_stacked_bar_20communities_proportions_ExpGroup"),
            colour = treat_col, 
            width = 6.5, 
            height = 6.5)
            #split = "multiple")  # 根据需要调整split参数
#domain
stacked_bar(fre_community20, X="agglomerate_20", select="new_domain", fill = "Tumour", OUTPUT_PATH, paste0(date, "_Rphenograph20_250526_stacked_bar_20communities_proportions_TumourDomain"),
            position = "count", colour = domain_col, width = 8, height = 6.5, split = "multiple")
# 


##check without self defined function , just -----------------

p = ggplot(Rphenograph20_250408, aes(
  x = reorder(community, ExpGroup == "WT", FUN = mean),
  fill = as.factor(ExpGroup)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  xlab("Community") +
  labs(fill = "") 

# 打印检查映射的中间图（检查基础ggplot是否正常）
print(p)

# 如果 position == "fill"，将显示比例
p_fill = p + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution")
print(p_fill)  # 打印显示比例的图

# # 如果 position == "count"，将显示计数
# p_count = p + geom_bar(stat = "count") + ylab("Cell count")
# print(p_count)  # 打印显示计数的图

# 先添加颜色映射并检查
p_coloured = p_fill + scale_fill_manual(values = treat_col)
print(p_coloured)  # 检查颜色是否正常

# 如果上一步正常，再添加 facet_wrap()
p_split = p_coloured + facet_wrap(~ExpGroup) + theme(strip.text = element_text(size = 16))
print(p_split)  # 检查拆分面板是否正常


# Commulative flow diagram--------------------------------------

neighber1<-Rphenograph20_250408
for (r in unique(neighber1$ROI_ID)) {
  cd = neighber1[neighber1$ROI_ID == r, ]
  
  tum_x_min = min(cd$Location_Center_X) 
  tum_x_max = max(cd$Location_Center_X) 
  tum_y_min = min(cd$Location_Center_Y) 
  tum_y_max = max(cd$Location_Center_Y) 
  
  tum_x_centre = (tum_x_min + tum_x_max)/2
  tum_y_centre = (tum_y_min + tum_y_max)/2
  
  neighber1[neighber1$ROI_ID == r, "scaled_X"] =
    (cd$Location_Center_X - tum_x_centre)/(tum_x_max - tum_x_centre)
  
  neighber1[neighber1$ROI_ID == r, "scaled_Y"] =
    (cd$Location_Center_Y - tum_y_centre)/(tum_y_max - tum_y_centre)
}

colnames(neighber1)

# Step 1: 筛选你要绘图的数据
df = neighber1[
  neighber1$scaled_Y <= 0.5 &
  neighber1$scaled_Y >= -0.5 &
  neighber1$domain == "TumDomain" &
  neighber1$ExpGroup == "WT", ]

# Step 2: 创建基础图层 + 灰色背景区域
p = ggplot(df, aes(x = scaled_X)) +
  annotate("rect", xmin = -Inf, xmax = -1.1, ymin = 0, ymax = Inf, alpha = 0.5, fill = "#F3F3F3") +
  annotate("rect", xmin = -1.1, xmax = -0.9, ymin = 0, ymax = Inf, alpha = 0.5, fill = "#979797") +
  annotate("rect", xmin = -0.9, xmax = 0.9,  ymin = 0, ymax = Inf, alpha = 0.5, fill = "#000000") +
  annotate("rect", xmin = 0.9, xmax = 1.1,   ymin = 0, ymax = Inf, alpha = 0.5, fill = "#979797") +
  annotate("rect", xmin = 1.1, xmax = Inf,   ymin = 0, ymax = Inf, alpha = 0.5, fill = "#F3F3F3")

# Step 3: 叠加 area 图
p = p + geom_area(aes(y = ..count.., fill = factor(community)), stat = "bin")
print(p)

# Step 4: 加标签、配色、美化
p = p +
  ylim(0, 4500) +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  guides(fill = guide_legend(title = "Community")) +
  xlab("Cross section through tissue") +
  ylab("Cell count") +
  scale_fill_manual(values = random_cols)
print(p)
# Step 5: 添加肿瘤方向的标注（推荐加上）
p = p +
  annotate("text", x = -1.3, y = 4000, label = "Tumour edge", size = 3, hjust = 0) +
  annotate("text", x = 0, y = 4000, label = "Tumour center", size = 3, hjust = 0.5) +
  annotate("text", x = 1.3, y = 4000, label = "Tumour edge", size = 3, hjust = 1) +
  annotate("segment", x = -1.2, xend = -1.0, y = 3800, yend = 3800,
           arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  annotate("segment", x = 1.2, xend = 1.0, y = 3800, yend = 3800,
           arrow = arrow(length = unit(0.2, "cm")), color = "black")
print(p)
# Step 6: 保存图像
ggsave(
  plot = p,
  device = "png",
  width = 7,
  height = 5.5,
  bg = 'white',
  dpi = 300,
  path = OUTPUT_PATH,
  filename = paste0(date, "_dataset1_cumulative_flow_20communities_Vehicle_onWTscale.png")
)


```


```{r}
# Fig 2b - MRTX

cumulative_flow(neighb1, "MRTX", 4.5, out_fig2, "_dataset1_cumulative_flow_18communities_MRTX_4500threshold.png")

# Supp Fig 2a
com18 = neighb1 %>% filter(agglom18_average == 18)
cumulative_flow(com18, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com18_Vehicle_onMRTX.png")
cumulative_flow(com18, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com18_MRTX_4500threshold.png")

# Supp Fig 2b
com10 = neighb1 %>% filter(agglom18_average == 10)
cumulative_flow(com10, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com10_Vehicle_onMRTX.png")
cumulative_flow(com10, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com10_MRTX_4500threshold.png")

# Supp Fig 2c
com2 = neighb1 %>% filter(agglom18_average == 2)
cumulative_flow(com2, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com2_Vehicle_onMRTX.png")
cumulative_flow(com2, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com2_MRTX_4500threshold.png")

# Supp Fig 2d
com3 = neighb1 %>% filter(agglom18_average == 3)
cumulative_flow(com3, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com3_Vehicle_onMRTX.png")
cumulative_flow(com3, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com3_MRTX_4500threshold.png")


#
```

# heatmap
```{r}
# Community proportions heatmap
community_prop_heatmap = function(data,
                                  ROI,
                                  name,
                                  out_path){
  
  #Generate table for heatmap that contains community size and proportion in treatment groups per ROI
  cd_prop = data.frame()
  for (r in unique(data[,ROI])){
    print(r)
    cd = data[which(data[,ROI] == r),]
    cd_prop[r, "ExpGroup"] = as.character(unique(cd$ExpGroup))
    cd_prop[r, ROI] = r
    total_cc = dim(cd)[1]
    for (cl in as.character(sort(unique(data$community)))){
      print(cl)
      cd_prop[r, cl] = dim(cd[which(cd$community == cl),])[1]/total_cc
    }
  }
  
  names(cd_prop)
  # Heatmap of cell type proportions
  heatmap_prop = cd_prop[,]
  row.names(heatmap_prop) = cd_prop[,ROI]
  head(heatmap_prop)
  heatmap_prop2 = heatmap_prop[,c(3:ncol(heatmap_prop))]
  heatmap_prop2 = as.matrix(heatmap_prop2)
  col1= colorRampPalette(brewer.pal(8, "Blues"))(25)
  col3 = list(ExpGroup = c("WT" = "#F8766D", "CCR2KO" = "#00BFC4"))
  dd = hclust(dist(heatmap_prop2), method = "average")
  ha = HeatmapAnnotation(ExpGroup = heatmap_prop$ExpGroup,
                         which = "row",
                         col = c(col3))
  

 # pdf(file=paste(out_path, date, name, sep = ""), width=9, height=5)
  ht = Heatmap(t(scale(t(heatmap_prop2))), name = "heatmap", col = col1, cluster_rows = dd, row_dend_side = "left", right_annotation = ha)
  print(ht) 
  #draw(ht)
 # dev.off()
}
community_prop_heatmap(fre_community20, "ROI_ID", "_neighb1_heatmap_community_proportions_per_ROI.png", OUTPUT_PATH)

```

# correlation Heatmap 
```{r}
# Create list of cell type pairs 
marker_list = function(celltypes){
  
  # Create df of each unique cell type 
  ct_pairs = data.frame('name1' = select)
  # Repeat each row in df by number of rows 
  ct_pairs = ct_pairs %>% slice(rep(1:n(), each = nrow(ct_pairs)))
  # Create another df on unique cell types and repeat block by number of rows
  pairs2 = data.frame('name2' = select)
  pairs2 = bind_rows(replicate(nrow(pairs2), pairs2, simplify = FALSE))
  # Merge columns into 1 df 
  ct_pairs = cbind(ct_pairs, pairs2)
  
  return(ct_pairs)
}

# Correlation of cell type pairs per community
community_correlations = function(data,
                                  name1,
                                  name2,
                                  df){
  
  # Create data frame with calculated proportions
  # of each cell type in each community
  prop = data %>%
    dplyr::select('source_cluster', 'community') %>%
    dplyr::group_by(community) %>%
    dplyr::mutate(community_count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(source_cluster %in% c(name1, name2)) %>%
    dplyr::group_by(community) %>%
    dplyr::summarise(community_count = unique(community_count),
                     count1 = sum(source_cluster == name1),
                     count2 = sum(source_cluster == name2),
                     prop1 = (count1/community_count)*100,
                     prop2 = (count2/community_count)*100
    )
  
  df$name1[row] = name1
  df$name2[row] = name2
  df$corr[row] = cor(prop$prop1, prop$prop2, method = "pearson")
  df$p_val[row] = cor.test(prop$prop1, prop$prop2, method = "pearson")$p.value
  
  return(df)
}

# Correlation of cell type pairs per ROI 
ct_correlations = function(data,
                           name1,
                           name2,
                           df){
  
  prop = data %>% 
    dplyr::select('source_cluster', 'ROI_ID') %>% 
    dplyr::group_by(ROI_ID) %>% 
    dplyr::mutate(roi_count = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(source_cluster %in% c(name1, name2)) %>% 
    dplyr::group_by(ROI_ID) %>% 
    dplyr::summarise(roi_count = unique(roi_count),
                     count1 = sum(source_cluster == name1),
                     count2 = sum(source_cluster == name2),
                     prop1 = (count1/roi_count)*100,
                     prop2 = (count2/roi_count)*100
    )
  
  df$name1[row] = name1
  df$name2[row] = name2
  df$corr[row] = cor(prop$prop1, prop$prop2, method = "pearson")
  df$p_val[row] = cor.test(prop$prop1, prop$prop2, method = "pearson")$p.value
  
  return(df)
}

select = c("Dendritic cells", "Dendritic cells CD103", "T cell CD4", "Tregs", "T cell CD8", "Neutrophils",
          "Macrophages", "B cells", "Fibroblasts", "EPTumour", "Tumour")

# Create a list of cell type pairs 
pairs = marker_list(select)

#### Figure 2e ####

# Create new df for saving cor value and p-value
corr_df = data.frame(matrix(nrow = nrow(pairs), ncol = 4))
names(corr_df) = c("name1", "name2", "corr", "p_val")

## Correlation of cell types per community ##
for(row in 1:nrow(pairs)){
  print(paste0(pairs$name1[row], ", ", pairs$name2[row]))
  corr_df = community_correlations(fre_community20, pairs$name1[row], pairs$name2[row], corr_df)
}

# Cluster cell types based on similar correlations ,sometimes doesnt work,use tge next code 
# corr_df_wide = corr_df %>% 
#   select(name1, name2, corr) %>% 
#   pivot_wider(names_from = name2, values_from = corr, values_fill = list(name1 = 0)) %>% 
#   column_to_rownames(var = "name1")
# 
# corr_df_wide = as.matrix(corr_df_wide)

#alternativecode based on tibber not dplyr
corr_df_wide <- corr_df %>% 
  select(name1, name2, corr) %>% 
  pivot_wider(
    names_from = name2, 
    values_from = corr, 
    values_fill = 0  # Fix: Use `0` instead of `list(name1 = 0)`
  ) %>% 
  tibble::column_to_rownames(var = "name1") %>%  # Explicitly use tibble's function
  as.matrix()


clust = hclust(dist(corr_df_wide), method = "average")
# Only plot the heatmap to see the order of the cell types following clustering 
heatmap.2(corr_df_wide, cluster_rows = clust,trace = "none", margins = c(10,10))
#dev.off()

# Order based on results from clust (plotted as heatmap for visualisation)
order = c("Dendritic cells", "Dendritic cells CD103", "T cell CD4", "Tregs", "T cell CD8", "Neutrophils",
          "Macrophages", "B cells", "Fibroblasts", "EPTumour", "Tumour")


# Add significance column to corr_df
corr_df$sig = ""
corr_df <- within(corr_df, sig[p_val < 0.05] <- '*')
corr_df <- within(corr_df, sig[p_val < 0.01] <- "**")
corr_df <- within(corr_df, sig[p_val > 0 & p_val < 0.001] <- "***")
corr_df <- within(corr_df, sig[corr == 0 & p_val == 0] <- "" )
# Replace NA values with zeros
corr_df[is.na(corr_df)] <- 0

# Create new DF for correlations with self 
corr_with_self = filter(corr_df, name1 == name2)
corr_df = filter(corr_df, !(name1 == name2))

# Remove repeated correlations : be careful ,not run now until clear know which is repeared correlations
# 添加一列标识成对信息（忽略顺序）
corr_df$pair_id = apply(corr_df[, c("name1", "name2")], 1, function(x) paste(sort(x), collapse = "_"))

# 去重（只保留每对组合的一次出现）
corr_df = corr_df[!duplicated(corr_df$pair_id), ]

# 删除 pair_id 列（可选）
corr_df$pair_id = NULL
#去掉自相关行
corr_df = corr_df[corr_df$name1 != corr_df$name2, ]
## this can be replaced
# corr_df = corr_df[-c(20,56,57,74,75,79,65,66,70,71,38,39,42:44,1,2,4,6:8,
#                       28:31,33:35,46:53,82:90),]
# Plot 
p = ggplot(corr_df,
           aes(x = factor(name1), y = factor(name2), fill = corr)) + 
  geom_tile()+
  theme_classic()+
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1)) + 
  geom_tile(data = corr_with_self, fill = "red") + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(label=sig), size = 5) + 
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0)) +
  scale_x_discrete(limits = order, NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0), limits=rev(order))
p
ggsave(plot = p, device = 'png', height = 5, width = 6.5, path = OUTPUT_PATH,
       filename = paste0(date, "_dataset1_cellType_correlation_matrix_perCommunity_half_reordered.png"))

#########################################################################################################
```
## Pearson correlation within each ROI ,not community . this need to combined with last chunk .
```{r}
#### Figure : Correlation of cell types per ROI ####

# 1. 相关性计算
roi_corr_df = data.frame(matrix(nrow = nrow(pairs), ncol = 4))
names(roi_corr_df) = c("name1", "name2", "corr", "p_val")

for(row in 1:nrow(pairs)){
  roi_corr_df = ct_correlations(fre_community20, pairs$name1[row], pairs$name2[row], roi_corr_df)
}

# 2. 标记显著性
roi_corr_df$sig = ""
roi_corr_df <- within(roi_corr_df, sig[p_val < 0.05] <- '*')
roi_corr_df <- within(roi_corr_df, sig[p_val < 0.01] <- "**")
roi_corr_df <- within(roi_corr_df, sig[p_val < 0.001] <- "***")
roi_corr_df <- within(roi_corr_df, sig[corr == 0 & p_val == 0] <- "" )
roi_corr_df[is.na(roi_corr_df)] <- 0

# 3. 去掉重复和自相关
roi_corr_df$pair_id = apply(roi_corr_df[, c("name1", "name2")], 1, function(x) paste(sort(x), collapse = "_"))
roi_corr_df = roi_corr_df[!duplicated(roi_corr_df$pair_id), ]
roi_corr_df$pair_id = NULL
roi_corr_with_self = filter(roi_corr_df, name1 == name2)
roi_corr_df = filter(roi_corr_df, !(name1 == name2))

# 4. 画图
p_roi = ggplot(roi_corr_df,
           aes(x = factor(name1), y = factor(name2), fill = corr)) + 
  geom_tile()+
  theme_classic()+
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1)) + 
  geom_tile(data = roi_corr_with_self, fill = "red") + 
  geom_text(aes(label=sig), size = 5) + 
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title = element_blank()) +
  scale_x_discrete(limits = order, NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0), limits=rev(order))

p_roi

# 5. 保存
ggsave(plot = p_roi, device = 'png', height = 5, width = 6.5, path = OUTPUT_PATH,
       filename = paste0(date, "_dataset1_cellType_correlation_matrix_perROI_half_reordered.png"))


#########################################################################################################
```




# map one community to one ROI-->check 2R4 ,maybe need to manually change name 
```{r}
#### Figure 3c ####

# Spatial location of communities on each ROI
######Did not include : "MOC2_CCR2KO_2R_4_ROI1_x1", "MOC2_CCR2KO_2R_4_ROI1_x2", because the code runing doesnt work. 
# Spatial location of communities on each ROI
# colours =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
       #      "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF",  "TCRgdt" = "#336666","Macrophages" = "#4E79A7FF",
#               "Macrophages F480" = "#4E79A7FF", "Macrophages CD68" = "#5F9EA0", "Macrophages CD11b" = "#0072B2", "Neutrophils" = "#A0CBE8FF",  "T cell CD4" = "#B07AA1FF","Monocytes"= "#40E0D0",  "T cell CD8" = "#FF9D9AFF", "Tregs" = "#CC6666", "Tumour" = "#704850FF", "Tumour EPCAM" = "#705","Tumour PancK" = "#752")


# Write CSV
#write.csv(celldata, paste0(OUTPUT_PATH, date, "celldata_vs.csv"), row.names = FALSE)

# Ensure the levels of the factor are set correctly
date<- format(Sys.Date(), "%y%m%d")


```
## 250408 Works
```{r}
# Prepare labels for plotting
images1 = c("MOC2_CCR2KO_1B_3_ROI1_1",  "MOC2_CCR2KO_1B_3_ROI2_-_s3", "MOC2_CCR2KO_1B_3_ROI2_2",  "MOC2_CCR2KO_1R_2_ROI1_1" ,  "MOC2_CCR2KO_NM_1_ROI1_1","MOC2_WT_1B_2_ROI1_1", "MOC2_WT_1R_4_ROI1_x1", "MOC2_WT_2R_3_ROI1_1",  "MOC2_WT_NM_1_ROI1_1","MOC2_CCR2KO_2R_4_ROI1")
out_dir = paste0(OUTPUT_PATH, "Community_map/")
dir.create(out_dir)
colnames(fre_community20)

# Table of comparable communities and labels - created in powerpoint 
fre_community20$vs_column<-"others"
fre_community20[which(fre_community20$community == '15'), "vs_column"] <- "15"
fre_community20[which(fre_community20$community == '1'), "vs_column"] <- "1"
fre_community20[which(fre_community20$community == '2'), "vs_column"] <- "2"
fre_community20[which(fre_community20$community == '7'), "vs_column"] <- "7"
fre_community20[which(fre_community20$community == '6'), "vs_column"] <- "6"
colnames(fre_community20)
unique(fre_community20$vs_column)

# Define variables

vs_column <- c("15", "1", 2, "7", "6")
vs_column_col <- c("15" = "#704850FF", "1" = "#FFF4A4FF", "2" = "#40E0D0", "7" = "green", "6" = "#FF9933", "Other" = "#E1E1E1")
fre_community20$vs_column <- factor(fre_community20$vs_column, levels = vs_column)
unique(fre_community20$vs_column)


nbclusters <- vs_column
labels <- as.data.frame(vs_column_col[1:5], stringsAsFactors = FALSE)
labels <- rownames_to_column(labels, "cluster")
colnames(labels)[2] <- 'colour'
add <- data.frame(cluster = c("101", "102"), colour = c("black", "white"), stringsAsFactors = FALSE)
labels <- rbind(add, labels)

# Define the cell_map function
cell_map <- function(ROI, files, data, clusters, path, labels) {
  
  # Filter the data for the specific ROI
  d <- data[data[[files]] == ROI, ]
  print(paste("ROI:", ROI, "Cells in ROI:", nrow(d)))
  # Construct file paths for TIFF images
  base_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Data/Spatial_input/"
  TIFFnm <- file.path(base_path, "all_cells_mask", paste0(ROI, "_all_cells_mask.tiff"))
  TIFFol <- file.path(base_path, "Cells_outline", paste0(ROI, "_Cells_outline.tiff"))
  
  # Print file paths for debugging
  print(paste("TIFFnm:", TIFFnm))
  print(paste("TIFFol:", TIFFol))
  
  # Read TIFF images
  if (!file.exists(TIFFnm)) {
    stop(paste("File not found:", TIFFnm))
  }
  if (!file.exists(TIFFol)) {
    stop(paste("File not found:", TIFFol))
  }
  
  TIFF <- readTIFF(TIFFnm)
  TIFF2 <- readTIFF(TIFFol)
  
  # Convert TIFF images to long format data frames
  to_long_df <- function(TIFF) {
    names(TIFF) <- seq_len(length(TIFF))
    TIFF <- melt(TIFF, id.vars = c(row(TIFF), names(TIFF)))
    names(TIFF)[1] <- "y"
    names(TIFF)[2] <- "x"
    TIFF
  }
  
  Outline <- TIFF2
  ROI1 <- TIFF
  ROI1 <- to_long_df(ROI1)
  Outline <- to_long_df(Outline)
  ROI1$unique_px_ID <- seq_len(nrow(ROI1))
  ROI2 <- ROI1
  n <- 2
  
  # Choose the plotting of neighbor clusters
for (cl in nbclusters) {
    print(cl)
    cluster_xy = d[which(d$vs_column == cl), c("Location_Center_X", "Location_Center_Y", "vs_column")]
    names(cluster_xy)[names(cluster_xy) == 'vs_column'] <- 'cluster' 
    # If cluster cl is not found in the image/ROI, remove it from the labels list
    if (dim(cluster_xy)[1] == 0) {
        print(paste("cluster ", cl, " is not found in ", ROI, sep = ""))
        next
    } else {
        cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
        cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
        names(cluster_xy) = c("x", "y", "cluster")
        colours_in_mask <- inner_join(ROI1, cluster_xy[, c("x", "y", "cluster")])
        
        # Check if 'value' in colours_in_mask is NULL or NA before processing
        if (is.null(colours_in_mask$value) | any(is.na(colours_in_mask$value))) {
            print("Warning: colours_in_mask$value contains NULL or NA values")
            next
        }
        
        min_val = min(unique(colours_in_mask$value))
        print(min_val)
        colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min_val)), ]
        ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
        ROI2$name[ROI2$value == n] = first(cluster_xy$cluster)
        n = n + 1
    }
}

print("now we come to the bit that is suspicious")

# Check if ROI2$name is NULL or contains NA values before recoding
if (is.null(ROI2$name)) {
    ROI2$name <- rep("NA", nrow(ROI2))  # Default to "NA" if NULL
}

# Handle any NA values in ROI2$name
ROI2$name[is.na(ROI2$name)] <- "NA"  # Replace NA with "NA"

# Recode ROI2$name using dplyr::recode
ROI2$name <- dplyr::recode(ROI2$name,
                           "1" = "15",
                           "2" = "1",
                           "3" = "2",
                           "4" = "7",
                           "5" = "6",
                           .default = "NA")

 # ("Tumour", "T cell CD8", "EPTumour", "Tregs", "Dendritic cells")
  background = ROI2[which(ROI2$value < 1),"unique_px_ID"]
  ROI2[which(Outline$value == 1 ),"value"] = 1
  ROI2[which(ROI2$unique_px_ID %in% background),"value"] = 0
  ROI2$name[ROI2$value == 0] = 101
  ROI2$name[ROI2$value == 1] = 102
  
  
  ## Subset labels based on communities present in that image
  labels = labels[which(labels$cluster %in% unique(ROI2$name)),]
  
  # Change 101 & 102 values to blank
  labels$cluster[labels$cluster == 101 | labels$cluster == 102] <- " "
  
  p = ggplot(ROI2,aes(x=x,y=-y, fill=as.factor(value))) +
    geom_raster() +
    theme_void() +
    theme(legend.title=element_blank(),
          legend.text = element_text(colour = "black")) +
  scale_fill_manual(values = alpha(labels$colour, 1), labels = labels$cluster) +
  ggtitle(basename(ROI))+
  theme(plot.title = element_text(hjust = 0.5, size = 10)) 
  ###################
  
  # filename = paste(ROI, ".pdf", sep = "")
  # ggsave(plot = p, device = "pdf", width=5.6, height=5, dpi=300, path = path, filename = filename)
  # print("plot saved")
   filename = paste(ROI, ".jpeg", sep = "")
  ggsave(plot = p, device = "jpeg", width=5.6, height=5, dpi=300, path = path, filename = filename)
  print("plot saved")
  return(p)
}

# Apply cell_map function to each ROI in images1
for (ROI in images1) {
  print(ROI)
  p <- cell_map(ROI, "filename", fre_community20, vs_column, out_dir, labels)
}
```

# step by step check the error for making map .   the code works--> go back to loop  
```{r}
## looking for error
# Step 1: 载入数据并设置路径
ROI <- "MOC2_CCR2KO_NM_1_ROI1_1"  # 输入你想测试的ROI名字
files <- "filename"   # 你的数据列名
data <- fre_community20   # 数据集
clusters <- vs_column   # 设置你要检查的簇
path <- out_dir   # 输出路径
labels <- labels   # 标签数据

# Step 2: Subset数据
d <- data[data[[files]] == ROI, ]
colnames(data)
colnames(d)
print(paste("ROI:", ROI, "Cells in ROI:", nrow(d)))

base_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Data/Figures_input/"
TIFFnm <- file.path(base_path, "all_cells_mask", paste0(ROI, "_all_cells_mask.tiff"))
TIFFol <- file.path(base_path, "Cells_outline", paste0(ROI, "_Cells_outline.tiff"))

# 检查文件是否存在
if (!file.exists(TIFFnm)) stop(paste("File not found:", TIFFnm))
if (!file.exists(TIFFol)) stop(paste("File not found:", TIFFol))

# Step 3: 读取TIFF文件
TIFF <- readTIFF(TIFFnm)
TIFF2 <- readTIFF(TIFFol)
print(paste("TIFF range:", range(TIFF)))
print(paste("Outline TIFF range:", range(TIFF2)))

# Step 4: 转换TIFF为长格式数据框
to_long_df <- function(TIFF) {
  names(TIFF) <- seq_len(length(TIFF))
  TIFF <- melt(TIFF, id.vars = c(row(TIFF), names(TIFF)))
  names(TIFF)[1] <- "y"
  names(TIFF)[2] <- "x"
  TIFF
}

Outline <- TIFF2
ROI1 <- TIFF
ROI1 <- to_long_df(ROI1)
Outline <- to_long_df(Outline)
ROI1$unique_px_ID <- seq_len(nrow(ROI1))
ROI2 <- ROI1
n <- 2

# Step 5: 映射簇
# Step 1: Loop over each cluster
for (cl in clusters) {
  print(paste("Processing cluster:", cl))

  # Step 2: Ensure the column exists, update with the correct column name if needed
  cluster_xy <- d[d$vs_column == cl, c("Location_Center_X", "Location_Center_Y", "vs_column")]
  
  # Step 3: Handle case where cluster is not found
  if (nrow(cluster_xy) == 0) {
    print(paste("Cluster", cl, "not found"))
    next
  }
  
  # Step 4: Rename and round the coordinates
  names(cluster_xy)[3] <- "cluster"  # due to i delete  object number from raw code ,so this should be 3 not 4
  cluster_xy$x <- round(cluster_xy$Location_Center_X)
  cluster_xy$y <- round(cluster_xy$Location_Center_Y)
  cluster_xy <- cluster_xy[, c("x", "y", "cluster")]

  # Step 5: Join with ROI1 and filter out minimum value
  colours_in_mask <- inner_join(ROI1, cluster_xy[, c("x", "y", "cluster")])
  print(paste("Matched pixels:", nrow(colours_in_mask)))

  # Step 6: If matched pixels are found, update the value and name
  if (nrow(colours_in_mask) > 0) {
    min_val <- min(unique(colours_in_mask$value))
    colours_in_mask <- colours_in_mask[colours_in_mask$value != min_val, ]
    ROI2$value[ROI1$value %in% colours_in_mask$value] <- n
    ROI2$name[ROI2$value == n] <- as.character(first(cluster_xy$cluster))
    n <- n + 1
  }
}


# Step 6: 重新编码簇名称
print("Step 4: Recode cluster names")
print("Unique before recode:")
print(unique(ROI2$name))
ROI2$name <- as.character(ROI2$name)
table(ROI2$name)
  ROI2$name <- dplyr::recode(ROI2$name,
                             "1" = "EP/Mac/T",
                             "2" = "Tum/Mac/T",
                             "3" = "EP/Mac/Neu/T",
                             "4" = "EP/T",
                             "5" = "DC/Mac/T",
                             .default = "NA")

background <- ROI2[ROI2$value < 1, "unique_px_ID"]
ROI2[Outline$value == 1, "value"] <- 1
ROI2[ROI2$unique_px_ID %in% background, "value"] <- 0
ROI2$name[ROI2$value == 0] <- 101
ROI2$name[ROI2$value == 1] <- 102

# Step 7: 准备标签
print("Step 5: Prepare labels")
print("Value table:")
print(table(ROI2$value))
print("Name table:")
print(table(ROI2$name))

labels <- labels[labels$cluster %in% unique(ROI2$name), ]
labels$cluster[labels$cluster %in% c(101, 102)] <- " "

# Step 8: 绘图
p <- ggplot(ROI2, aes(x = x, y = -y, fill = as.factor(value))) +
  geom_raster() +
  theme_void() +
  theme(legend.title = element_blank(), legend.text = element_text(color = "black")) +
  scale_fill_manual(values = alpha(labels$colour, 1), labels = labels$cluster) +
  ggtitle(basename(ROI)) +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

filename <- paste0(ROI, "_debug.jpeg")
ggsave(plot = p, device = "jpeg", width = 5.6, height = 5, dpi = 300, path = path, filename = filename)
print(paste("Plot saved to", file.path(path, filename)))

# Apply cell_map function to each ROI in images1
for (ROI in images1) {
  print(ROI)
  p <- cell_map(ROI, "filename", fre_community20, vs_column2, out_dir, labels)
}

```

#Stichting
```{r}
gc()  # 触发垃圾回收，释放内存
library(magick)
# 定义 stitch_images 函数
stitch_images <- function(folder_path, n_rows = 2, n_cols = 5, image_format = "jpeg", output_format = "jpeg") {
  
  jpeg_files <- list.files(folder_path, pattern = paste0("\\.", image_format, "$"), full.names = TRUE)
  
  images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
  num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图
  
  for (i in 1:num_batches) {
    start_idx <- (i - 1) * images_per_combined + 1
    end_idx <- min(i * images_per_combined, length(jpeg_files))
    
    batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
    
    images <- lapply(batch_files, image_read)  # 读取当前批次的图片
    
    combined_image <- image_montage(
      image_join(images),
      tile = paste0(n_cols, "x", n_rows),  # 拼接的行列数
      geometry = "x500",  # 每个图片的大小
      bg = "white"  # 背景色
    )
    
    output_path <- file.path(folder_path, paste0("combined_image_", i, ".", output_format))
    image_write(combined_image, path = output_path, format = output_format)
    cat("Image saved at:", output_path, "\n")
  }
}

# 使用该函数时可以自行设置 n_rows 和 n_cols
stitch_images("/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Spatial_output/Community/Community_map", n_rows = 2, n_cols = 5,image_format = "jpeg")


```

#### Functions ####

# Plotting expression of markers across communities & treatment groups - box plot 

```{r}
## merge celldata 
data1$top5<-"others"
data1[which(data1$community == '15'), "top5"] <- "15"
data1[which(data1$community == '1'), "top5"] <- "1"
data1[which(data1$community == '2'), "top5"] <- "2"
data1[which(data1$community == '7'), "top5"] <- "7"
data1[which(data1$community == '6'), "top5"] <- "6"

colnames(data1)
unique(data1$top5)
top5_col <- c("15" = "#704850FF", "1" = "#FFF4A4FF", "2" = "#40E0D0", "7" = "green", "6" = "#FF9933", "Other" = "#E1E1E1")

#### Functions ####
# Plotting expression of markers across communities & treatment groups - box plot 
expression_boxplot = function(data, select_ct, X, marker, dir, name, strip_text){
  
  p = ggplot(transform(data[which(data$top5 %in% c("15", "1", "2", "7", "6") & data$annotation2 %in% select_ct),],
                       ExpGroup = factor(ExpGroup, levels = c("WT", 'CCR2KO'))),
             aes(x = ExpGroup, y = get(X), fill = top5))  + 
    geom_boxplot() + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          strip.text = element_text(size=strip_text, face = "bold"),
          legend.text = element_text(size = 16),
          legend.title = element_blank(),
          legend.position = "right") +
    ylab(paste0(marker, " expression")) +
    scale_y_continuous(trans='log2', labels = label_number(accuracy = 0.0001)) +
    scale_fill_manual(values= top5_col[1:5]) +
    facet_wrap(annotation2~.)
  p
  ggsave(plot = p, device = "png", width=8.5, height=6, dpi=300, path = dir,
         filename = paste0(date, name))
  print(p)
}

expression_boxplot(data1, c("Dendritic cells", "Dendritic cells CD103"), "MI_PDL1", "PD-L1", OUTPUT_PATH, "_data1_top5_PDL1_expression_DCs.png", 12.8)
expression_boxplot(data1, c("Dendritic cells", "Dendritic cells CD103"), "MI_CD86", "CD86", OUTPUT_PATH, "_data1_top5_CD86_expression_DCs.png", 12.8)
expression_boxplot(data1, c("Dendritic cells", "Dendritic cells CD103"), "MI_CD80", "CD80", OUTPUT_PATH, "_data1_top5_CD80_expression_DCs.png", 12.8)
expression_boxplot(data1, c("Dendritic cells", "Dendritic cells CD103"), "MI_CD24", "CD24", OUTPUT_PATH, "_data1_top5_CD24_expression_DCs.png", 12.8)
expression_boxplot(data1, c("Dendritic cells", "Dendritic cells CD103"), "MI_MHCcII", "MHCII", OUTPUT_PATH, "_data1_top5_MHCII_expression_DCs.png", 12.8)



expression_boxplot(data1, c("Macrophages"), "MI_CD86", "CD86", OUTPUT_PATH, "_data1_top5_CD86_expression_MACs.png", 12.8)
expression_boxplot(data1, c("Macrophages"), "MI_CD206", "CD206", OUTPUT_PATH, "_data1_top5_CD20_expression_MACs.png", 12.8)
expression_boxplot(data1, c("Macrophages"), "MI_CD68", "CD68", OUTPUT_PATH, "_data1_top5_CD68_expression_MACs.png", 12.8)
expression_boxplot(data1, c("Macrophages"), "MI_CD11b", "CD11b", OUTPUT_PATH, "_data1_top5_CD11b_expression_MACs.png", 12.8)
expression_boxplot(data1, c("Macrophages"), "MI_MHCcII", "MHCII", OUTPUT_PATH, "_data1_top5_MHCII_expression_MACs.png", 12.8)
expression_boxplot(data1, c("Macrophages"), "MI_PDL1", "PD-L1", OUTPUT_PATH, "_data1_top5_PD-L1_expression_MACs.png", 12.8)
expression_boxplot(data1, c("Macrophages"), "MI_F480", "F480", OUTPUT_PATH, "_data1_top5_PD-L1_expression_MACs.png", 12.8)

expression_boxplot(data1, c("Tumour","EPTumour"), "MI_CD44", "CD44", OUTPUT_PATH, "_data1_top5_CD44_expression_Tums.png", 12.8)
expression_boxplot(data1, c("Tumour","EPTumour"), "MI_CD24", "CD24", OUTPUT_PATH, "_data1_top5_CD24_expression_Tums.png", 12.8)
expression_boxplot(data1, c("Tumour","EPTumour"), "MI_PDL1", "PD-L1", OUTPUT_PATH, "_data1_top5_PD-L1_expression_Tums.png", 12.8)
expression_boxplot(data1, c("Tumour","EPTumour"), "MI_Ki67", "Ki67", OUTPUT_PATH, "_data1_top5_Ki67_expression_Tums.png", 12.8)
expression_boxplot(data1, c("Tumour","EPTumour"), "MI_casp3", "Caspase3", OUTPUT_PATH, "_data1_top5_Casp3_expression_Tums.png", 12.8)
expression_boxplot(data1, c("Tumour","EPTumour"), "MI_pS6", "pS6", OUTPUT_PATH, "_data1_top5_pS6_expression_Tums.png", 12.8)

stitch_images("/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Spatial_output/Community/Community_boxplot/tumour", n_rows = 2, n_cols = 3,image_format = "png")
gc() 
stitch_images("/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Spatial_output/Community/Community_boxplot/dc", n_rows = 2, n_cols = 3,image_format = "png")
gc() 
stitch_images("/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Spatial_output/Community/Community_boxplot/mac", n_rows = 2, n_cols = 4,image_format = "png")
```

# Plotting expression of markers across communities & treatment groups - density plot 
```{r}
# Plotting expression of markers across communities & treatment groups - density plot 
expression_density = function(data, select_ct, X, marker, xlow, xhigh, dir, name, w, h){

  p = ggplot(transform(data[which(data$top5 %in% c("15", "1", "2", "7", "6")  & data$annotation2 %in% select_ct
                                     & data[,X] != '-Inf'),],
                       ExpGroup = factor(ExpGroup, levels = c("WT", 'CCR2KO'))),
             aes(x = get(X), y = top5, fill = ExpGroup))  + 
    geom_density_ridges(scale = 0.9, alpha = 0.8, quantile_lines = T, quantile_fun = median) + 
    scale_y_discrete(expand = c(0.05, 0)) + 
    xlim(xlow,xhigh) +  
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size=20, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          legend.position = "right") +
    ylab("Community") +
    scale_fill_manual(values= treat_col) 
  if(length(select_ct) == 1){
    p = p +  xlab(paste0(marker, " expression on ", select_ct)) 
  } else (
    p = p + xlab(paste0(marker, " expression")) 
  )
  p  
  ggsave(plot = p, device = "png", dpi=300, width=w, height=h, path = dir,
         filename = paste0(date, name))
  print(p)
}
data1$log_PD1 = log2(data1$MI_PD1)
expression_density(data1, "T cell CD8", "log_PD1", "PD-1", -7, 3, OUTPUT_PATH, "_dataset1_top5_CD8_PD1_expression_treatment_comparison_density.png", 5.5, 6)
data1$log_TIM3 = log2(data1$MI_TIM3)
expression_density(data1, "T cell CD8", "log_TIM3", "TIM3", -7, 3, OUTPUT_PATH, "_dataset1_top5_CD8_TIM3_expression_treatment_comparison_density.png", 5.5, 6)
data1$log_LAG3 = log2(data1$MI_LAG3)
expression_density(data1, "T cell CD8", "log_LAG3", "LAG3", -7, 3, OUTPUT_PATH, "_dataset1_top5_CD8_LAG3_expression_treatment_comparison_density.png", 5.5, 6)
# expression_density(data1, c("Dendritic cells", "Dendritic cells CD103", "Macrophages type 1", "Macrophages type 2"), "log_CXCL9", "CXCL9", -10, 5, out_fig4,
#                    "_dataset2_top5_CXCL9_expression_treatment_comparison_density.png", 6, 7)
stitch_images("/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Spatial_output/Community/Community_Density", n_rows = 2, n_cols = 4,image_format = "png")
```
# Plotting expression of markers to better understand tumour phenotype

```{r}
# Plotting expression of markers to better understand tumour phenotype
tum_expr = function(data, Y, com, cols, marker, dir, name){

  p = ggplot(data[which(data$cellType == 'Tumour' & data$top5 %in% com & data$treatment == "MRTX"),],
             aes(x = factor(top5), y = get(Y), colour = factor(top5))) + 
    geom_boxplot(size = 2) + 
    theme_classic() + 
    scale_color_manual(values = cols) +
    scale_y_continuous(trans='log2', labels = label_number(accuracy = 0.0001)) + 
  xlab("") + 
    ylab(paste0("Tumour cell ", marker, " expression")) + 
    labs(color = "") + 
    theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 22),
          axis.text.x = element_blank(),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 22))
  p
  ggsave(plot = p, device = "png", width=5, height=5, dpi=300, path = dir,
         filename = paste0(date, name))
}

#########################################################################################################


################
#### Fig 4a ####



#### Fig 4c ####

# Calculated logged values of CXCL9 expression 
neighb2$log_CXCL9 = log2(neighb2$MI_CXCL9)
expression_density(neighb2, c("Dendritic cells", "Dendritic cells CD103", "Macrophages type 1", "Macrophages type 2"), "log_CXCL9", "CXCL9", -10, 5, out_fig4,
                   "_dataset2_top5_CXCL9_expression_treatment_comparison_density.png", 6, 7)

################
#### Fig 4d ####

# May need to change this to look at all 5 communities 
tum_expr(neighb1, "MI_Ki67", c("T/DC", "T/M2_1", "T/M2_2"), top5_col[3:5], "Ki67", out_fig4, "_dataset1_top5_tumour_ki67_expression_MRTXonly.png")



#### Fig 4e ####

tum_expr(neighb1, "MI_casp3", c("T/DC", "T/M2_1", "T/M2_2"), top5_col[3:5], "C-casp3", out_fig4, "_dataset1_top5_tumour_casp3_expression_MRTXonly.png")


#### Fig 4f ####

# Calculated logged values of PD-1 expression 
neighb1$log_PD1 = log2(neighb1$MI_PD1)
neighb2$log_PD1 = log2(neighb2$MI_PD1)

expression_density(neighb1, "T cells CD8", "log_PD1", "PD-1", -7, 3, out_fig4, "_dataset1_top5_CD8_PD1_expression_treatment_comparison_density.png", 5.5, 6)


#####################
#### Supp Fig 4k ####

## Percentage of PD-1+ CD8 T cells that are LAG3+ 
p = neighb2 %>%
  dplyr::select(cellType, top5, MI_PD1, MI_LAG3, treatment) %>%
  dplyr::filter(
    cellType == "T cells CD8" &
      MI_PD1 >= 0.5 &
      top5 %in% c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2") & treatment == "MRTX"
  ) %>%
  dplyr::group_by(top5) %>%
  dplyr::mutate(pd1_cd8_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(top5) %>%
  dplyr::summarise(
    pd1_cd8_count = unique(pd1_cd8_count),
    lag3_pd1_cd8_count = sum(MI_LAG3 > 0.5),
    percent_pos = (lag3_pd1_cd8_count / pd1_cd8_count) * 100
  ) %>%
  ggplot(aes(
    x = top5,
    y = percent_pos,
    fill = top5
  )) +
  geom_col(colour = "Black", position = "dodge") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    axis.text.x = element_text(size = 22, angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  scale_colour_manual(values = treat_col) +
  xlab("") +
  ylab("PD-1+ CD8 T cells: percentage LAG3+") +
  ylim(0, 40)
p
ggsave(
  plot = p,
  device = "png",
  width = 4.5,
  height = 7,
  dpi = 300,
  path = out_supp4,
  filename = paste(
    date,
    "_dataset2_PD1CD8_percetnage_LAG3pos_community_treatment.png",
    sep = ""
  )
)


```




