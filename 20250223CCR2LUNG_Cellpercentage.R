#Cell percentage

library(ggpubr) 

BASE = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG"

# Name of input and output directories:
INPUT_PATH = file.path(BASE, "Data/Annotation_input")
OUTPUT_PATH = file.path(BASE, "Results/Cellperc_output")
celldata = read.csv(paste0(INPUT_PATH, "/250526_celldata_clean.csv"))

celldata$LAG3exp <- ifelse(celldata$MI_LAG3 > 1, "LAG3+", "LAG3-") # alittle bit more, difficult, some signal seems come from Ki67, but keep the Lag3 signal which coexpress with T cell
celldata$Ki67exp <- ifelse(celldata$MI_Ki67 >= 0.2, "Ki67+", "Ki67-") #find a balance
celldata$Casp3exp <-ifelse(celldata$MI_Casp3>=0.7,"Casp3+","Casp3-")
celldata$PDL1exp <- ifelse(celldata$MI_PDL1 >= 0.3, "PDL1+", "PDL1-")
celldata$pS6exp <- ifelse(celldata$MI_pS6 >= 0.3, "pS6+", "pS6-") 
celldata$MHCIIexp <- ifelse(celldata$MI_MHCII >= 0.25, "MHCII+", "MHCII-") 
celldata$PD1exp <- ifelse(celldata$MI_PD1 >= 0.5, "PD1+", "PD1-") 
celldata$CD86exp <- ifelse(celldata$MI_CD86 >= 0.8, "CD86+", "CD86-") #set a little bit strict , consider noise
celldata$CD206exp <- ifelse(celldata$MI_CD206>= 0.3, "CD206+", "CD206-") 
celldata$PVRexp <- ifelse(celldata$MI_PVR >= 0.6, "PVR+", "PVR-") 
celldata$TIM3exp <- ifelse(celldata$MI_TIM3 >= 0.8, "TIM3+", "TIM3-") 
celldata$Thy1exp<- ifelse(celldata$MI_Thy1_2 >= 0.2, "Thy1+", "Thy1-") 
celldata$CXCL9exp<- ifelse(celldata$MI_CXCL9 >= 0.5, "CXCL9+", "CXCL9-") 

# Calculate the percentage of cell types 

str(celldata)
celldata$ExpGroup <- as.factor(celldata$ExpGroup)
celldata$ROI_ID <- as.factor(celldata$ROI_ID)
celldata$MouseID <- as.factor(celldata$MouseID)
# str(celldata)
# # Remove the unnamed column
# celldata <- celldata[, !colnames(celldata) == "X"]
# # Verify the column removal
# colnames(celldata)

###################################################
# what is include here:
##1 all cell number   2. T ,Mac ,Mon, Neu cells in immune cells  3.CD4,8 Treg percentage in T cell population. 
#activation marker and exaustion marker. 


# stat_250523_ROI = celldata %>%
#   group_by(ExpGroup,MouseID,ROI_ID) %>%
  
stat_250523_Mice = celldata %>%
   group_by(ExpGroup,MouseID) %>%
  dplyr::summarise(
    totaln = n(),
    Tumourn = sum(annotation=="Tumour"),
    Tcelln = sum(annotation %in% c("T cell CD4","T cell CD8","Tregs")),
    CD4n = sum(annotation=="T cell CD4"),
    CD8n = sum(annotation=="T cell CD8"),
    Tregn= sum(annotation=="Tregs"),
    Bcelln= sum(annotation=='B cells'),
    Neutn = sum(annotation=="Neutrophils"),
    AMn = sum(annotation=="Alvelor Macrophages"),
    IMn = sum(annotation=="Interstitial Macrophages"),
    Maccelln = sum(annotation %in% c("Alvelor Macrophages","Interstitial Macrophages")),
    Monon=sum(annotation=='Monocytes'),
    Fibn = sum(annotation=="Fibroblasts"),
    Endon=sum(annotation=="Endothelium"),
    Epin=sum(annotation=="Epithelium"),
    cDC1n = sum(annotation =="Dendritic cells CD103"),
    DCothern = sum(annotation =="Dendritic cells other"),
    DCn =sum (annotation %in% c("Dendritic cells CD103","Dendritic cells other")),
    immunecellsn = sum(annotation %in% c("T cell CD4","T cell CD8","Tregs","Neutrophils","Alvelor Macrophages","Interstitial Macrophages","B cells","Dendritic cells CD103","Dendritic cells other",'Monocytes')),
    
    # T cell surface marker : PD1  TIM3  LAG3 Ki67
    PD1Tcelln = sum(PD1exp=="PD1+" & annotation %in% c("T cell CD4","T cell CD8","Tregs")),
    PD1CD4Tcelln = sum(PD1exp=="PD1+" & annotation=="T cell CD4"),
    PD1CD8Tcelln = sum(PD1exp=="PD1+" & annotation=="T cell CD8"),
    PD1Tregn = sum(PD1exp=="PD1+" & annotation=="Tregs"),
    PD1Neutn = sum(PD1exp=="PD1+" & annotation=="Neutrophils"),
    PD1Monon = sum(PD1exp=="PD1+" & annotation=="Monocytes"),
    
    LAG3Tcelln = sum(LAG3exp=="LAG3+" & annotation %in% c("T cell CD4","T cell CD8","Tregs")),
    LAG3CD4Tcelln = sum(LAG3exp=="LAG3+" & annotation=="T cell CD4"),
    LAG3CD8Tcelln = sum(LAG3exp=="LAG3+" & annotation=="T cell CD8"),
    LAG3Tregn = sum(LAG3exp=="LAG3+" & annotation=="Tregs"),
    
    
    TIM3Tcelln = sum(TIM3exp=="TIM3+" & annotation %in% c("T cell CD4","T cell CD8","Tregs")),
    TIM3CD4Tcelln = sum(TIM3exp=="TIM3+" & annotation=="T cell CD4"),
    TIM3CD8Tcelln = sum(TIM3exp=="TIM3+" & annotation=="T cell CD8"),
    TIM3Tregn = sum(TIM3exp=="TIM3+" & annotation=="Tregs"),
    
    
    Ki67Tcelln = sum(Ki67exp=="Ki67+" & annotation %in% c("T cell CD4","T cell CD8","Tregs")),
    Ki67CD4Tcelln = sum(Ki67exp=="Ki67+" & annotation=="T cell CD4"),
    Ki67CD8Tcelln = sum(Ki67exp=="Ki67+" & annotation=="T cell CD8"),
    Ki67Tregn = sum(Ki67exp=="Ki67+" & annotation=="Tregs"),
    Ki67Endon= sum(Ki67exp=="Ki67+" & annotation=="Endothelium"),
    # CD80 CD86 on Neutrophils , PD1 PDL1 TIM3 LAG3 Ki67 MHCII
    
    #MHCII, CD86, CD80  PD1 PDL1 MHCII PVR on myeloid cells
    CD86Maccelln = sum(CD86exp=="CD86+" & annotation %in% c("Alvelor Macrophages","Interstitial Macrophages")),
    CD206Maccelln = sum(CD206exp=="CD206+" & annotation %in% c("Alvelor Macrophages","Interstitial Macrophages")),
    CD86AMn = sum(CD86exp=="CD86+" & annotation == "Alvelor Macrophages"),
    CD86IMn = sum(CD86exp=="CD86+" & annotation == "Interstitial Macrophages"),
    CD86cDC1n = sum(CD86exp=="CD86+" & annotation=="Dendritic cells CD103"),
    CD86DCothern = sum(CD86exp=="CD86+" & annotation=="Dendritic cells other"),
    CD86DCn = sum(CD86exp=="CD86+" & annotation %in% c("Dendritic cells CD103","Dendritic cells other")),
    CD86Monon = sum(CD86exp=="CD86+" & annotation=="Monocytes"),
    CD86Neutn = sum(CD86exp=="CD86+" & annotation=="Neutrophils"),
    
    pS6Maccelln = sum(pS6exp=="pS6+" & annotation %in% c("Alvelor Macrophages","Interstitial Macrophages")),
    pS6AMn = sum(pS6exp=="pS6+" & annotation == "Alvelor Macrophages"),
    pS6IMn = sum(pS6exp=="pS6+" & annotation == "Interstitial Macrophages"),
    pS6Monon = sum(pS6exp=="pS6+" & annotation=="Monocytes"),
    pS6Neutn = sum(pS6exp=="pS6+" & annotation=="Neutrophils"),
    pS6Tumourn= sum(pS6exp=="pS6+" & annotation=="Tumour"),
    
    MHCIIMaccelln = sum(MHCIIexp=="MHCII+" & annotation %in% c("Alvelor Macrophages","Interstitial Macrophages")),
    MHCIIAMn = sum(MHCIIexp=="MHCII+" & annotation == "Alvelor Macrophages"),
    MHCIIIMn = sum(MHCIIexp=="MHCII+" & annotation == "Interstitial Macrophages"),
    MHCIIcDC1n = sum(MHCIIexp=="MHCII+" & annotation=="Dendritic cells CD103"),
    MHCIIDCothern = sum(MHCIIexp=="MHCII+" & annotation=="Dendritic cells other"),
    MHCIIDCn = sum(MHCIIexp=="MHCII+" & annotation %in% c("Dendritic cells CD103","Dendritic cells other")),
    
    PDL1Maccelln = sum(PDL1exp=="PDL1+" & annotation %in% c("Alvelor Macrophages","Interstitial Macrophages")),
    PDL1AMn = sum(PDL1exp=="PDL1+" & annotation == "Alvelor Macrophages"),
    PDL1IMn = sum(PDL1exp=="PDL1+" & annotation == "Interstitial Macrophages"),
    PDL1Neutn = sum(PDL1exp=="PDL1+" & annotation=="Neutrophils"),
    PDL1cDC1n = sum(PDL1exp=="PDL1+" & annotation=="Dendritic cells CD103"),
    PDL1DCothern = sum(PDL1exp=="PDL1+" & annotation=="Dendritic cells other"),
    PDL1DCn = sum(PDL1exp=="PDL1+" & annotation %in% c("Dendritic cells CD103","Dendritic cells other")),
    PDL1Monon = sum(PDL1exp=="PDL1+" & annotation=="Monocytes"),
    PDL1Tumourn = sum(PDL1exp=="PDL1+" & annotation=="Tumour"),
    PDL1Endon = sum(PDL1exp=="PDL1+" & annotation=="Endothelium"),
    
    PVRMaccelln = sum(PVRexp=="PVR+" & annotation %in% c("Alvelor Macrophages","Interstitial Macrophages")),
    PVRAMn = sum(PVRexp=="PVR+" & annotation == "Alvelor Macrophages"),
    PVRIMn = sum(PVRexp=="PVR+" & annotation == "Interstitial Macrophages"),
    PVRMonon = sum(PVRexp=="PVR+" & annotation=="Monocytes"),
    PVRTumourn = sum(PVRexp=="PVR+" & annotation=="Tumour"),
    PVRFibn = sum(PVRexp=="PVR+" & annotation=="Fibroblasts"),
    #  casp3  Ki67  PVR  pS6 on tumour
    
    Casp3Tumourn = sum(Casp3exp=="Casp3+" & annotation=="Tumour"),
    Ki67Tumourn = sum(Ki67exp=="Ki67+" & annotation=="Tumour")
  )%>%
  group_by(ExpGroup,MouseID) %>%
 # group_by(ExpGroup,MouseID,ROI_ID) %>%
  mutate(          
    Tcelln_immunecellsn = Tcelln/immunecellsn,
    Bcelln_immunecells = Bcelln/immunecellsn,
    Neutn_immunecells = Neutn/immunecellsn,
    Maccelln_immunecells= Maccelln/immunecellsn, 
    AMn_immunecells = AMn/immunecellsn,
    IMn_immunecells = IMn/immunecellsn,
    cDC1n_immunecells= cDC1n/immunecellsn,
    DCothern_immunecells= DCothern/immunecellsn,
    Monon_immunecells= Monon/immunecellsn, 
    Fibn_totaln =Fibn/totaln,
    Endon_totaln = Endon/totaln,
    Epin_totaln =Epin/totaln,
    
    # Each sub cell type %  in each cell type 
    CD4n_Tcelln = CD4n/Tcelln,
    CD8n_Tcelln = CD8n/Tcelln,
    Tregn_Tcelln = Tregn/Tcelln,
    CD4n_immunecellsn = CD4n/immunecellsn,
    CD8n_immunecellsn = CD8n/immunecellsn,
    Tregn_immunecellsn = Tregn/immunecellsn,
    
    cDC1n_immunecellsn = cDC1n/immunecellsn,
    DCothern_immunecellsn = DCothern/immunecellsn,
    Maccelln_immunecellsn = Maccelln/immunecellsn,
    Monon_immunecellsn = Monon/immunecellsn,
    Neutn_immunecellsn = Neutn/immunecellsn,
    
    # Marker expression relative to cell type
    PD1Tcelln_Tcelln = PD1Tcelln/Tcelln,
    PD1CD8Tcelln_CD8n =  PD1CD8Tcelln/CD8n,
    PD1CD4Tcelln_CD4n = PD1CD4Tcelln/CD4n,
    PD1Tregn_Tregn = PD1Tregn/Tregn,
    cDC1n_DCn = cDC1n/DCn,
    DCothern_DCn = DCothern/DCn,
    
    CD86cDC1n_cDC1n = CD86cDC1n/cDC1n,
    CD86DCothern_DCothern = CD86DCothern/DCothern,
    CD86DCn_DCn =CD86DCn/DCn,
    MHCIIcDC1n_cDC1n = MHCIIcDC1n/cDC1n,
    MHCIIDCothern_DCothern = MHCIIDCothern/DCothern,
    MHCIIDCn_DCn =MHCIIDCn/DCn,
    PDL1cDC1n_cDC1n = PDL1cDC1n/cDC1n,
    PDL1DCothern_DCothern = PDL1DCothern/DCothern,
    PDL1DCn_DCn =PDL1DCn/DCn,
    
    
    MHCIIMaccelln_Maccelln =MHCIIMaccelln/Maccelln,
    CD86Maccelln_Maccelln =CD86Maccelln/Maccelln,
    CD206Maccelln_Maccelln =CD206Maccelln/Maccelln,
    AMn_Maccelln=AMn/Maccelln,
    IMn_Maccelln=IMn/Maccelln,
    
    
    PDL1Maccelln_Maccelln = PDL1Maccelln/Maccelln,
    PVRMaccelln_Maccelln =PVRMaccelln/Maccelln,
    pS6Maccelln_Maccelln =pS6Maccelln/Maccelln,
    
    Casp3Tumourn_Tumourn = Casp3Tumourn/Tumourn,
    Ki67Tumourn_Tumourn =  Ki67Tumourn/Tumourn,
    PDL1Tumourn_Tumourn = PDL1Tumourn/Tumourn,
    PVRTumourn_Tumourn =  PVRTumourn/Tumourn,
    pS6Tumourn_Tumourn =  pS6Tumourn/Tumourn,
    
    
    PD1Neutn_Neutn =PD1Neutn/Neutn,
    PDL1Neutn_Neutn =PDL1Neutn/Neutn,
    CD86Neutn_Neutn =CD86Neutn/Neutn,
    
    
    PD1Monon_Monon =PD1Monon/Monon,
    PDL1Monon_Monon =PDL1Monon/Monon,
    CD86Monon_Monon =CD86Monon/Monon,
    pS6Monon_Monon =pS6Monon/Monon,
  ) %>%
  ungroup()


plot.boxplot.ROI = function(celldata, comparisons = list(c("WT", "CCR2KO"))) {
  numeric_cols <- celldata %>% select(where(is.numeric)) %>% colnames()
  for (factor in numeric_cols) {
    factor_sym <- ensym(factor)  
    is_percentage <- max(celldata[[factor]], na.rm = TRUE) <= 1.5
    p <- ggplot(celldata, aes(x = ExpGroup, y = !!factor_sym, colour = MouseID, fill = MouseID)) +
      #geom_boxplot(alpha = 0.5) +  # 透明度 0.5 避免遮挡
      geom_dotplot(method = "histodot", binaxis = 'y', stackdir = 'center', dotsize = 0.8) +
      geom_text(aes(label = paste0(ROI_ID,MouseID)), position = position_jitter(width = 0.2, height = 0), 
                size = 2, hjust = -0.3) +
      (if (is_percentage) scale_y_continuous(labels = scales::percent) else scale_y_continuous()) + 
      scale_x_discrete(limits = c("WT", "CCR2KO")) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                         label = "p.signif", size = 6, fontface = "bold") +
      theme_classic() +
      ggtitle(factor) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.background = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.text.x = element_text(size = 17, face = "bold", colour = "black"),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size = 5))
    
    date <- Sys.Date()
    ggsave(filename = paste0(date, "_Boxplot_", factor, ".jpeg"),
           path = file.path(OUTPUT_PATH, "Boxplots_percentage/ROI"),
           plot = p, bg = "white",
           width = 2100, height = 1700, units = "px")
  }
}

my_comparisons <- list(c("WT", "CCR2KO"))
plot.boxplot.ROI(celldata = stat_250523_ROI, comparisons = my_comparisons)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)
plot.boxplot.Mice <- function(celldata, comparisons = list(c("WT", "CCR2KO"))) {

  
  numeric_cols <- celldata %>% select(where(is.numeric)) %>% colnames()
  
  # 初始化结果表格
  results <- data.frame(
    Marker = character(),
    p_value = numeric(),
    Method = character(),
    P_less_0.05 = logical(),
    P_less_0.01 = logical(),
    stringsAsFactors = FALSE
  )
  
  for (factor in numeric_cols) {
    # 判断是否是百分比列（值不超过1.5）
    if (max(celldata[[factor]], na.rm = TRUE) > 1.5) next
    
    factor_sym <- rlang::ensym(factor)
    
    # Wilcoxon 检验
    group1 <- celldata %>% filter(ExpGroup == comparisons[[1]][1]) %>% pull(!!factor_sym)
    group2 <- celldata %>% filter(ExpGroup == comparisons[[1]][2]) %>% pull(!!factor_sym)
    
    if (length(group1) > 0 && length(group2) > 0 && length(unique(c(group1, group2))) > 1) {
      test_result <- wilcox.test(group1, group2)
      p_val <- test_result$p.value
      p05 <- p_val < 0.05
      p01 <- p_val < 0.01
      
      results <- rbind(results, data.frame(
        Marker = factor,
        p_value = p_val,
        Method = "wilcox.test",
        P_less_0.05 = p05,
        P_less_0.01 = p01
      ))
    } else {
      next  # 若数据不足或恒定，跳过
    }
    
    # 绘图部分
    p <- ggplot(celldata, aes(x = ExpGroup, y = !!factor_sym, colour = MouseID, fill = MouseID)) +
      geom_dotplot(method = "histodot", binaxis = 'y', stackdir = 'center', dotsize = 0.8) +
      geom_text(aes(label = MouseID), position = position_jitter(width = 0.2, height = 0), 
                size = 2, hjust = -0.3) +
      scale_x_discrete(limits = c("WT", "CCR2KO")) +
      (if (max(celldata[[factor]], na.rm = TRUE) <= 1.5)
        scale_y_continuous(labels = percent)
       else scale_y_continuous()) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                         label = "p.signif", size = 6, fontface = "bold") +
      theme_classic() +
      ggtitle(factor) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.background = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.text.x = element_text(size = 17, face = "bold", colour = "black"),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size = 5))
    
    date <- Sys.Date()
    ggsave(filename = paste0(date, "_Boxplot_", factor, ".jpeg"),
           path = file.path(OUTPUT_PATH, "Boxplots_percentage/Mice"),
           plot = p, bg = "white",
           width = 2100, height = 1700, units = "px")
  }
  
  # 保存表格结果
  write.csv(results, file.path(OUTPUT_PATH, "Boxplots_percentage", "Wilcoxon_Results.csv"), row.names = FALSE)
  
  # 打印显著 marker
  cat("显著 marker (p < 0.05):\n")
  print(results$Marker[results$P_less_0.05])
}
my_comparisons <- list(c("WT", "CCR2KO"))
plot.boxplot.Mice(celldata = stat_250523_Mice, comparisons = my_comparisons)


# plot.boxplot.Mice = function(celldata, comparisons = list(c("WT", "CCR2KO"))) {
#   numeric_cols <- celldata %>% select(where(is.numeric)) %>% colnames()
#   for (factor in numeric_cols) {
#     factor_sym <- ensym(factor)  
#     is_percentage <- max(celldata[[factor]], na.rm = TRUE) <= 1.5
#     p <- ggplot(celldata, aes(x = ExpGroup, y = !!factor_sym, colour = MouseID, fill = MouseID)) +
#       #geom_boxplot(alpha = 0.5) +  # 透明度 0.5 避免遮挡
#       geom_dotplot(method = "histodot", binaxis = 'y', stackdir = 'center', dotsize = 0.8) +
#       geom_text(aes(label = paste0(MouseID)), position = position_jitter(width = 0.2, height = 0), 
#                 size = 2, hjust = -0.3) +
#       (if (is_percentage) scale_y_continuous(labels = scales::percent) else scale_y_continuous()) + 
#       scale_x_discrete(limits = c("WT", "CCR2KO")) +
#       stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
#                          label = "p.signif", size = 6, fontface = "bold") +
#       theme_classic() +
#       ggtitle(factor) +
#       theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#             plot.background = element_blank(),
#             axis.title = element_blank(),
#             axis.text.y = element_text(size = 14, colour = "black"),
#             axis.text.x = element_text(size = 17, face = "bold", colour = "black"),
#             legend.position = "none",
#             legend.title = element_blank(),
#             legend.text = element_text(size = 5))
#     
#     date <- Sys.Date()
#     ggsave(filename = paste0(date, "_Boxplot_", factor, ".jpeg"),
#            path = file.path(OUTPUT_PATH, "Boxplots_percentage/Mice"),
#            plot = p, bg = "white",
#            width = 2100, height = 1700, units = "px")
#   }
# }
# my_comparisons <- list(c("WT", "CCR2KO"))
# plot.boxplot.Mice(celldata = stat_250523_Mice, comparisons = my_comparisons)



#_____________________________________________________________________________________________________________________
## stackbarplot 

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
## set cell percentage based on domain
celldata <- celldata %>%
  mutate(across(c(domain, domain2), as.factor))
celldata_CCR2<-celldata[celldata$ExpGroup=="CCR2KO",]
celldata_WT<-celldata[celldata$ExpGroup=="WT",]
levels_annotation = c("T cell CD4", "T cell CD8", "Tregs", "B cells","Neutrophils", "Alvelor Macrophages", "Interstitial Macrophages","Monocytes", "Dendritic cells CD103",  "Dendritic cells other", "Endothelium", "Epithelium", "Fibroblasts","Tumour","Uncertain")

plot.celltypeperc <- function(celldata, title) {
  date <- str_sub(gsub("-", "", Sys.Date()), -6, -1)
  directory <- create.resultsdir(path = OUTPUT_PATH, name = "Stacked_barplots")
  filename <- paste0(title, ".jpeg")
  
  # Ensure domain levels are correctly set
  celldata$domain <- factor(celldata$domain, levels = c("Normal", "Interface", "Tumour", "n/a"))
  
  p <- ggplot(
    data = celldata[which(celldata$domain != "n/a"), ],
    aes(
      x = factor(domain, levels = c("Normal","Interface","Tumour")),
      fill = annotation  
    )
  ) +
    geom_bar(position = "fill", colour = "black") +
    scale_fill_manual(
      # values = c("T cell CD4" = "#B07AA1FF", "T cell CD8" = "#FF9D9AFF", "Tregs" = "#CC6666","B cells" = "#945931",  "Macrophages" = "#4E79A7FF", "Dendritic cells CD103" = "#FFCC66","Dendritic cells" = "#F26484", "Neutrophils" = "#A0CBE8FF","Endothelium" = "#FFEA42",  "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF","Tumour" = "#704850FF","Monocytes"="#336666")
      values =  c("B cells" = "#945931",  "Dendritic cells other" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
                  "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF",  
                  "Alvelor Macrophages" = "#4E79A7FF", "Interstitial Macrophages" = "#5F9EA0",  "Neutrophils" = "#A0CBE8FF",  "T cell CD4" = "#B07AA1FF","Monocytes"= "#40E0D0",  "T cell CD8" = "#FF9D9AFF", "Tregs" = "#CC6666", "Tumour" = "#704850FF", "Uncertain" = "#8C8C8C")
      
      
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 60, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 6),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.title = element_blank()
    )
  
  filename <- paste0(date, "_", title, ".jpeg")
  ggsave(
    device = "jpeg", width = 2000, height = 2000, units = "px", dpi = 300,
    path = directory, filename = filename
  )
  
  return(p)
}

# Call the function
plot.celltypeperc(celldata = celldata_CCR2, title = "Stacked_CCR2_Domain")
plot.celltypeperc(celldata = celldata_WT, title = "Stacked_WT_Domain")
# plot.celltypeperc(celldata = celldata[!celldata$annotation %in% c("Tumour", "Uncertain"),],
#                   title = "Stacked_without_tumour")
plot.celltypeperc(celldata = celldata_WT[!celldata_WT$annotation %in% c("Tumour", "Endothelium"),],
                  title = "Stacked_immunecells_WT_Domain")
plot.celltypeperc(celldata = celldata_CCR2[!celldata_CCR2$annotation %in% c("Tumour", "Endothelium"),],
                  title = "Stacked_immunecells_CCR2_Domain")


#x domain  y  annotation
plot.celltypeperc <- function(celldata, title) {
  date <- str_sub(gsub("-", "", Sys.Date()), -6, -1)
  directory <- create.resultsdir(path = OUTPUT_PATH, name = "Stacked_barplots")
  filename <- paste0(title, ".jpeg")
  
  # Ensure domain levels are correctly set
  celldata$domain <- factor(celldata$domain, levels = c("Normal", "Interface", "Tumour", "n/a"))
  
  p <- ggplot(
    data = celldata[which(celldata$domain != "n/a"), ],
    aes(
      x = factor(annotation, levels = levels_annotation),
      fill = domain  # Use domain directly, already properly leveled
    )
  ) +
    geom_bar(position = "fill", colour = "black") +
    scale_fill_manual(
      values = c(
        "Normal" = "#B07AA1FF",
        "Interface" = "#A0CBE8FF",
        "Tumour" = "#ABDDA4FF",
        "n/a" = "gray"
      )
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 60, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 6),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.title = element_blank()
    )
  
  filename <- paste0(date, "_", title, ".jpeg")
  ggsave(
    device = "jpeg", width = 2000, height = 2000, units = "px", dpi = 300,
    path = directory, filename = filename
  )
  
  return(p)
}



# Call the function
plot.celltypeperc(celldata = celldata_CCR2, title = "Stacked_domain_CCR2")
plot.celltypeperc(celldata = celldata_WT, title = "Stacked_domain_WT")
plot.celltypeperc(celldata = celldata_WT[!celldata_WT$annotation %in% c("Tumour", "Endothelium"),],
                  title = "Stacked_immunecells_WT_Domain")
plot.celltypeperc(celldata = celldata_CCR2[!celldata_CCR2$annotation %in% c("Tumour", "Endothelium"),],
                  title = "Stacked_immunecells_CCR2_Domain")



plot.celltypeperc <- function(celldata, title) {
  date <- str_sub(gsub("-", "", Sys.Date()), -6, -1)
  directory <- create.resultsdir(path = OUTPUT_PATH, name = "Stacked_barplots")
  filename <- paste0(title, ".jpeg")
  
  # Ensure domain levels are correctly set
  
  
  p <- ggplot(
    data = celldata[which(celldata$domain != "n/a"), ],
    aes(
      x = factor(ROI_ID, levels = unique(celldata$ROI_ID)),
      fill = annotation  
    )
  ) +
    geom_bar(position = "fill", colour = "black") +
    scale_fill_manual(
      values =  c("B cells" = "#945931",  "Dendritic cells other" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
                  "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF",  
                  "Alvelor Macrophages" = "#4E79A7FF", "Interstitial Macrophages" = "#5F9EA0",  "Neutrophils" = "#A0CBE8FF",  "T cell CD4" = "#B07AA1FF","Monocytes"= "#40E0D0",  "T cell CD8" = "#FF9D9AFF", "Tregs" = "#CC6666", "Tumour" = "#704850FF", "Uncertain" = "#8C8C8C")
      
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 60, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 6),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.title = element_blank()
    )
  
  filename <- paste0(date, "_", title, ".jpeg")
  ggsave(
    device = "jpeg", width = 2000, height = 2000, units = "px", dpi = 300,
    path = directory, filename = filename
  )
  
  return(p)
}

plot.celltypeperc(celldata = celldata, title = "Stacked_ROI")
plot.celltypeperc(celldata = celldata[!celldata$annotation %in% c("Tumour", "Endothelium","Epithelium"),],
                  title = "Stacked_immunecells_ROI")

## save the file.
date<- format(Sys.Date(), "%y%m%d")
output_file = paste0(OUTPUT_PATH,"/", date, "stat_250523_Mice.csv")#OUTPUT_PATH
output_file
write.csv(stat_250523_Mice, file = output_file, row.names = F) 
output_file = paste0(OUTPUT_PATH,"/", date, "stat_250523_ROI.csv")#OUTPUT_PATH
output_file
write.csv(stat_250523_ROI, file = output_file, row.names = F) 
#  linear mix model for ROI statitics ----------------------------------------------------------------------------------

#install.packages("lmerTest")
library(lmerTest)
library(lme4)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sjPlot)  # 可视化边际效应
library(purrr)
#install.packages("sjPlot")
names(stat_250523_ROI)

library(lme4)
library(lmerTest)

# 保证因子格式
stat_250523_ROI$ExpGroup <- factor(stat_250523_ROI$ExpGroup)

# 母类细胞字段定义
parent_terms <- c("totaln", "Tumourn", "Tcelln", "CD4n", "CD8n", "Tregn",
                  "Bcelln", "Neutn", "AMn", "IMn", "Maccelln", "Monon", "Fibn",
                  "Endon", "Epin", "cDC1n", "DCothern", "DCn", "immunecellsn")

# 识别所有比例列
all_cols <- names(stat_250523_ROI)
percent_cols <- all_cols[
  sapply(all_cols, function(x) {
    parts <- strsplit(x, "_")[[1]]
    length(parts) == 2 && parts[2] %in% parent_terms
  })
]

# 初始化结果表
results <- data.frame(
  Marker = character(),
  Estimate = numeric(),
  p_value = numeric(),
  UsedModel = character(),
  P_less_0.05 = logical(),
  P_less_0.01 = logical(),
  stringsAsFactors = FALSE
)

# 初始化跳过信息记录表
skipped <- data.frame(
  Marker = character(),
  Reason = character(),
  stringsAsFactors = FALSE
)

# 循环分析
for (col in percent_cols) {
  data_sub <- stat_250523_ROI[!is.na(stat_250523_ROI[[col]]), ]
  
  # 跳过值恒定列
  if (length(unique(data_sub[[col]])) <= 1) {
    skipped <- rbind(skipped, data.frame(Marker = col, Reason = "Only one unique value"))
    next
  }
  
  # 跳过只有一个组别的列
  if (length(unique(data_sub$ExpGroup)) < 2) {
    skipped <- rbind(skipped, data.frame(Marker = col, Reason = "Only one ExpGroup level after NA removal"))
    next
  }
  
  # 判断是否能使用 lmer
  mouse_counts <- table(data_sub$MouseID)
  has_repeated_mouse <- any(mouse_counts > 1)
  formula <- as.formula(paste0(col, " ~ ExpGroup"))
  
  tryCatch({
    if (has_repeated_mouse) {
      formula_lmer <- as.formula(paste0(col, " ~ ExpGroup + (1 | MouseID)"))
      model <- lmer(formula_lmer, data = data_sub)
      coef_summary <- summary(model)$coefficients
      method <- "lmer"
    } else {
      model <- lm(formula, data = data_sub)
      coef_summary <- summary(model)$coefficients
      method <- "lm"
    }
    
    group_row <- rownames(coef_summary)[grepl("^ExpGroup", rownames(coef_summary))]
    
    if (length(group_row) == 1) {
      estimate <- coef_summary[group_row, "Estimate"]
      p_val <- coef_summary[group_row, "Pr(>|t|)"]
      p05 <- p_val < 0.05
      p01 <- p_val < 0.01
      
      results <- rbind(results, data.frame(
        Marker = col,
        Estimate = estimate,
        p_value = p_val,
        UsedModel = method,
        P_less_0.05 = p05,
        P_less_0.01 = p01
      ))
    }
  }, error = function(e) {
    skipped <- rbind(skipped, data.frame(Marker = col, Reason = "Model fitting error"))
  })
}

# 保存结果和跳过信息
# 保存结果为CSV
write.csv(results, "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/Boxplots_percentage/250526LMM_pvalue_ROI_wholetissue.csv",row.names = FALSE)
write.csv(skipped, "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/Boxplots_percentage/LMM_skipped_columns.csv", row.names = FALSE)

# 输出显著marker
cat("显著 marker (p < 0.05):\n")
print(results$Marker[results$P_less_0.05])

cat("\n已跳过的 marker 数量:", nrow(skipped), "\n")



#  visualize the cell population _____________________________________________________________________________________________________________



##  New stackbarplot + statistics  based on dataset-------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(readr)
library(tidyr)
library(forcats)

##
# Filter out rows with domain == "N/A"
celldata <- celldata %>% filter(!is.na(new_domain), new_domain != "n/a", !is.na(annotation))
celldata$new_domain <- factor(celldata$new_domain, levels = c("Tumour", "Interface", "Normal"))
celldata$ExpGroup <- factor(celldata$ExpGroup, levels = c("WT", "CCR2KO"))
annotation_order <- c("T cell CD4", "T cell CD8", "Tregs", "B cells", 
                      "Dendritic cells CD103", "Dendritic cells other", 
                      "Monocytes", "Alvelor Macrophages", "Interstitial Macrophages",  "Neutrophils", 
                      "Fibroblasts", "Endothelium", "Epithelium", "Tumour", "Uncertain")

# Mouse-level Analysis
#-----------------------------#
# Mouse-level Analysis
#-----------------------------#
mouse_data <- celldata %>%
  group_by(MouseID, ExpGroup, new_domain, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(MouseID, ExpGroup, new_domain) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plot (mouse-level)
mouse_plot <- mouse_data %>%
  group_by(MouseID, ExpGroup, annotation, new_domain) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%
  group_by(ExpGroup, annotation) %>%
  mutate(total = sum(prop)) %>%
  group_by(ExpGroup, annotation) %>%
  mutate(prop = prop / total)  # Normalize to 1 within each annotation
mouse_plot$annotation <- factor(mouse_plot$annotation, levels = annotation_order)

if (nrow(mouse_plot) > 0 && length(unique(mouse_plot$ExpGroup)) > 0) {
  p_mouse <- ggplot(mouse_plot, aes(x = annotation, y = prop, fill = new_domain)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~ExpGroup) +
    theme_minimal() +
    ylab("Proportion") +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  print(p_mouse)
  ggsave("/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/StackedBarPlots/MouseLevel/stacked_bar_mouse_level.pdf", plot = p_mouse, width = 10, height = 6)
} else {
  message("No data available for mouse_plot.")
}

# Wilcoxon Test (mouse-level)
mouse_stats <- data.frame()

for (ann in unique(mouse_data$annotation)) {
  for (dom in unique(mouse_data$new_domain)) {
    sub <- mouse_data %>% filter(annotation == ann, new_domain == dom)
    if (length(unique(sub$ExpGroup)) == 2) {
      test <- wilcox.test(prop ~ ExpGroup, data = sub)
      mouse_stats <- rbind(mouse_stats, data.frame(
        annotation = ann,
        domain = dom,
        method = "wilcox.test",
        p_value = test$p.value,
        significant = ifelse(test$p.value < 0.05, TRUE, FALSE)
      ))
    }
  }
}

write_csv(mouse_stats, "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/StackedBarPlots/MouseLevel/stats_mouse_level.csv")

# Print significant results
cat("\n[ROI-Level] Significant results (p < 0.05):\n")
print(mouse_stats %>% filter(significant == TRUE) %>% select(annotation, domain))


# Plot dotplots for significant results,mouselevel ,wilvox-------
sig_mouse_stats <- mouse_stats %>% filter(significant == TRUE)

for (i in seq_len(nrow(sig_mouse_stats))) {
  ann <- sig_mouse_stats$annotation[i]
  dom <- sig_mouse_stats$domain[i]
  sub <- mouse_data %>% filter(annotation == ann, new_domain == dom)
  p <- ggplot(sub, aes(x = ExpGroup, y = prop, color = ExpGroup)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, fill = NA) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    theme_classic() +
    ylab("Proportion") +
    ggtitle(paste0("Mouse-Level: ", ann, " in ", dom)) +
    theme(axis.title.x = element_blank())
  
  ggsave(filename = paste0("/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/StackedBarPlots/MouseLevel/dotplot_", gsub(" ", "_", ann), "_", dom, ".pdf"),
         plot = p, width = 5, height = 4)
}



# ROI-level Analysis-------
#-----------------------------#
#-----------------------------#
roi_data <- celldata %>%
  group_by(ROI_ID, MouseID, ExpGroup, new_domain, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ROI_ID, new_domain) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Normalize within each annotation and ExpGroup for fill plot
roi_plot <- roi_data %>%
  group_by(ROI_ID, ExpGroup, annotation, new_domain) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%
  group_by(ExpGroup, annotation) %>%
  mutate(total = sum(prop)) %>%
  group_by(ExpGroup, annotation) %>%
  mutate(prop = prop / total)

if (nrow(roi_plot) > 0 && length(unique(roi_plot$ExpGroup)) > 0) {
  p_roi <- ggplot(roi_plot, aes(x = annotation, y = prop, fill = new_domain)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~ExpGroup) +
    theme_minimal() +
    ylab("Proportion") +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  print(p_roi)
  ggsave("/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/StackedBarPlots/ROILevel/stacked_bar_roi_level.pdf", width = 10, height = 6)
} else {
  message("No data available for roi_plot.")
}


# Linear Mixed Model (ROI-level)-----LMM ,seems not work yet --------
lmm_stats <- data.frame()

library(lme4)

library(lme4)
library(dplyr)
library(openxlsx)

library(lme4)
library(dplyr)

lmm_stats <- data.frame()
skipped_stats <- data.frame()

for (ann in unique(roi_data$annotation)) {
  for (dom in unique(roi_data$new_domain)) {
    sub <- roi_data %>% filter(annotation == ann, new_domain == dom)
    
    # 检查是否有两个组别
    if (length(unique(sub$ExpGroup)) < 2) {
      skipped_stats <- rbind(skipped_stats, data.frame(
        annotation = ann,
        domain = dom,
        reason = "Only one ExpGroup"
      ))
      next
    }
    
    # 尝试拟合模型
    tryCatch({
      model <- lmer(prop ~ ExpGroup + (1 | MouseID), data = sub)
      
      if (isSingular(model, tol = 1e-4)) {
        skipped_stats <- rbind(skipped_stats, data.frame(
          annotation = ann,
          domain = dom,
          reason = "Singular fit"
        ))
      } else {
        sm <- summary(model)
        row <- grep("^ExpGroup", rownames(sm$coefficients), value = TRUE)
        if (length(row) > 0) {
          p_val <- sm$coefficients[row, "Pr(>|t|)"]
          lmm_stats <- rbind(lmm_stats, data.frame(
            annotation = ann,
            domain = dom,
            method = "lmer",
            p_value = p_val,
            significant = p_val < 0.05
          ))
        } else {
          skipped_stats <- rbind(skipped_stats, data.frame(
            annotation = ann,
            domain = dom,
            reason = "No valid ExpGroup coefficient"
          ))
        }
      }
    }, error = function(e) {
      skipped_stats <- rbind(skipped_stats, data.frame(
        annotation = ann,
        domain = dom,
        reason = paste("Model error:", e$message)
      ))
    })
  }
}

# 保存为 CSV 文件
write.csv(lmm_stats, "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/StackedBarPlots/ROILevel/lmm_results.csv", row.names = FALSE)
write.csv(skipped_stats, "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Cellperc_output/StackedBarPlots/ROILevel/lmm_skipped.csv", row.names = FALSE)


# Print significant results
cat("\n[ROI-Level] Significant results (p < 0.05):\n")
print(lmm_stats %>% filter(significant == TRUE) %>% select(annotation, domain))
table(sub$ExpGroup)
table(sub$MouseID)
print(paste("Skipping:", ann, dom, "Rows:", nrow(sub)))
print(table(sub$ExpGroup))
print(table(sub$MouseID))




#MI--------------------------------------
#Baed on cell level ---------------------------------------------------------------
OUTPUT_PATH = file.path(BASE, "Results/Cellperc_output/MI_celltype/")
library(ggplot2)
library(scales)
library(ggpubr)

expgroup_colors <- c("WT" = "#FFF4A4FF", "CCR2KO" = "#40E0D0")

# 判断使用哪种检验
select_test_method <- function(df) {
  group_counts <- table(df$ExpGroup)
  method <- "wilcox.test"  # 默认
  reasons <- ""
  
  # 如果任一组样本数量 < 3，不做正态性检验
  if (any(group_counts < 3)) {
    reasons <- "Sample size < 3 in one or both groups: fallback to Wilcoxon"
    return(list(method = method, reason = reasons))
  }
  
  # 正态性检验（每组样本数量 >= 3）
  groups <- split(df[[2]], df$ExpGroup)
  try_results <- tryCatch({
    normal_results <- lapply(groups, function(g) {
      if (length(g) >= 3 && length(g) <= 5000) {
        shapiro.test(g)$p.value > 0.05
      } else {
        FALSE  # 超出范围也标记为非正态
      }
    })
    
    if (all(unlist(normal_results))) {
      method <- "t.test"
      reasons <- "Both groups pass Shapiro-Wilk normality test (p > 0.05)"
    } else {
      reasons <- "One or both groups fail normality test"
    }
    
    list(method = method, reason = reasons)
  }, error = function(e) {
    list(method = "wilcox.test", reason = paste("Shapiro test failed:", e$message))
  })
  
  return(try_results)
}


# 主函数
expression_boxplot <- function(celldata, select_ct, X, marker, dir, name, strip_text, image_format = "tiff") {
  date_str <- format(Sys.Date(), "%Y%m%d")
  plotdata <- celldata[celldata$annotation %in% select_ct, ]
  plotdata$ExpGroup <- factor(plotdata$ExpGroup, levels = c("WT", "CCR2KO"))
  
  stats_table <- data.frame(
    annotation = character(),
    test_method = character(),
    reason = character(),
    p_value = numeric(),
    signif = character(),
    stringsAsFactors = FALSE
  )
  
  comparisons <- list(c("WT", "CCR2KO"))
  
  p <- ggplot(plotdata, aes(x = ExpGroup, y = get(X), fill = ExpGroup)) +
    geom_boxplot() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = strip_text, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "right"
    ) +
    ylab(paste0(marker, " expression")) +
    scale_y_continuous(trans = 'log2', labels = label_number(accuracy = 0.0001)) +
    scale_fill_manual(values = expgroup_colors) +
    facet_wrap(annotation ~ .)
  
  annots <- unique(plotdata$annotation)
  for (a in annots) {
    sub <- plotdata[plotdata$annotation == a, c("ExpGroup", X)]
    colnames(sub)[2] <- "value"
    sub <- na.omit(sub)
    
    method_info <- select_test_method(sub)
    pval <- tryCatch({
      if (method_info$method == "t.test") {
        t.test(value ~ ExpGroup, data = sub)$p.value
      } else {
        wilcox.test(value ~ ExpGroup, data = sub)$p.value
      }
    }, error = function(e) NA)
    
    signif_label <- if (is.na(pval)) "NA"
    else if (pval <= 0.0001) "****"
    else if (pval <= 0.001) "***"
    else if (pval <= 0.01) "**"
    else if (pval <= 0.05) "*"
    else "ns"
    
    stats_table <- rbind(stats_table, data.frame(
      annotation = a,
      test_method = method_info$method,
      reason = method_info$reason,
      p_value = pval,
      signif = signif_label
    ))
  }
  
  write.csv(stats_table, file = file.path(dir, paste0(date_str, name, "_stats.csv")), row.names = FALSE)
  
  # 自动计算每个 annotation 的最大值
  y_max_by_group <- aggregate(get(X) ~ annotation, data = plotdata, max)
  names(y_max_by_group) <- c("annotation", "y_pos")
  y_max_by_group$y_pos <- y_max_by_group$y_pos * 1.5  # 星号略高于最大值
  
  # 添加星号标注：每个子图分别添加
  for (i in 1:nrow(y_max_by_group)) {
    this_annot <- y_max_by_group$annotation[i]
    this_y <- y_max_by_group$y_pos[i]
    
    # p <- p + stat_compare_means(
    #   data = subset(plotdata, annotation == this_annot),
    #   aes(label = ..p.signif..),
    #   method = NULL,
    #   label = "p.signif",
    #   comparisons = list(c("Vehicle", "Necroptosis")),
    #   size = 5,
    #   label.y = this_y
    # )
  }
  
  
  ggsave(plot = p,
         filename = paste0(date_str, name, ".", image_format),
         path = dir,
         device = image_format,
         width = 8.5,
         height = 6,
         dpi = 300)
  
  print(p)
}

##if i have multiple marker
batch_plot_markers <- function(celldata, select_ct, marker_list, dir, strip_text = 12.8, image_format = "tiff") {
  for (marker in marker_list) {
    X <- paste0("MI_", marker)
    name <- paste0("_", marker, "_expression_Mac")
    expression_boxplot(
      celldata = celldata,
      select_ct = select_ct,
      X = X,
      marker = marker,
      dir = dir,
      name = name,
      strip_text = strip_text,
      image_format = image_format
    )
  }
}
## call 
marker_list <- c("CD86", "CD68", "Sirpa", "MHCII", "PDL1")

selected_cts <- c( "Alvelor Macrophages","Interstitial Macrophages","Monocytes",'Dendritic cells other','Dendritic cells CD103')

batch_plot_markers(celldata, selected_cts, marker_list, OUTPUT_PATH, image_format = "tiff")

######Based on MOuse level-------------------------------------------------
library(ggplot2)
library(ggpubr)
library(scales)

# 自定义颜色
expgroup_colors <- c("WT" = "#f6e58d", "CCR2KO" = "#7ed6df")

library(ggplot2)
library(dplyr)
library(ggpubr)
library(scales)

# Custom colors for experimental groups
expgroup_colors <- c("WT" = "#f6e58d", "CCR2KO" = "#7ed6df")

# Function for plotting and statistical testing for a single marker
expression_boxplot <- function(celldata, select_ct, X, marker, dir, name, strip_text = 12.8, image_format = "png") {
  date_str <- format(Sys.Date(), "%Y%m%d")
  
  # Step 1: Aggregate by mean per mouse per cell type
  agg_formula <- as.formula(paste(X, "~ MouseID + annotation + ExpGroup"))
  aggdata <- aggregate(agg_formula, data = celldata, FUN = mean, na.rm = TRUE)
  
  # Step 2: Filter selected cell types
  plotdata <- aggdata[aggdata$annotation %in% select_ct, ]
  plotdata$ExpGroup <- factor(plotdata$ExpGroup, levels = c("WT", "CCR2KO"))
  
  # Step 3: Compute y-axis position for p-value stars
  y_max_by_group <- aggregate(get(X) ~ annotation, data = plotdata, max)
  names(y_max_by_group) <- c("annotation", "y_pos")
  y_max_by_group$y_pos <- pmin(y_max_by_group$y_pos * 2.5, 100)
  
  # Step 4: Perform Wilcoxon tests
  stats_table <- data.frame(
    annotation = character(),
    test_method = character(),
    p_value = numeric(),
    signif = character(),
    stringsAsFactors = FALSE
  )
  compare_layers <- list()
  
  for (i in seq_len(nrow(y_max_by_group))) {
    this_annot <- y_max_by_group$annotation[i]
    subdata <- subset(plotdata, annotation == this_annot)
    
    if (length(unique(subdata$ExpGroup)) < 2) next
    
    test_res <- tryCatch({
      pval <- wilcox.test(get(X) ~ ExpGroup, data = subdata)$p.value
      method <- "wilcox.test"
      signif <- if (pval <= 0.0001) "****" else if (pval <= 0.001) "***"
      else if (pval <= 0.01) "**" else if (pval <= 0.05) "*" else "ns"
      
      stats_table <- rbind(stats_table, data.frame(
        annotation = this_annot,
        test_method = method,
        p_value = pval,
        signif = signif
      ))
      
      stat_compare_means(
        data = subdata,
        aes(label = ..p.signif..),
        method = method,
        label = "p.signif",
        comparisons = list(c("WT", "CCR2KO")),
        size = 5,
        label.y = y_max_by_group$y_pos[i]
      )
    }, error = function(e) NULL)
    
    if (!is.null(test_res)) {
      compare_layers[[length(compare_layers) + 1]] <- test_res
    }
  }
  
  # Step 5: Plot
  p <- ggplot(plotdata, aes(x = ExpGroup, y = get(X), fill = ExpGroup)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    facet_wrap(annotation ~ .) +
    scale_y_continuous(trans = 'log2', labels = label_number(accuracy = 0.0001)) +
    scale_fill_manual(values = expgroup_colors) +
    ylab(paste0(marker, " expression (mean per mouse)")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      strip.text = element_text(size = strip_text, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  if (length(compare_layers) > 0) {
    for (layer in compare_layers) {
      p <- p + layer
    }
  }
  
  # Step 6: Save plot
  ggsave(
    plot = p,
    filename = paste0(date_str, name, ".", image_format),
    path = dir,
    device = image_format,
    width = 8.5,
    height = 6,
    dpi = 300
  )
  
  # Step 7: Save stats table
  stats_path <- file.path(dir, paste0(date_str, name, "_stats.csv"))
  write.csv(stats_table, file = stats_path, row.names = FALSE)
  
  # Return stats for later use
  stats_table$marker <- marker
  return(stats_table)
}

# Batch version for multiple markers
batch_plot_markers <- function(celldata, select_ct, marker_list, dir, strip_text = 12.8, image_format = "tiff") {
  all_stats <- list()
  
  for (marker in marker_list) {
    X <- paste0("MI_", marker)
    name <- paste0("_", marker, "_mouse_level_expression")
    
    stats_table <- expression_boxplot(
      celldata = celldata,
      select_ct = select_ct,
      X = X,
      marker = marker,
      dir = dir,
      name = name,
      strip_text = strip_text,
      image_format = image_format
    )
    
    all_stats[[marker]] <- stats_table
  }
  
  all_stats_df <- do.call(rbind, all_stats)
  return(all_stats_df)
}

# Example input
marker_list <- c("CD86", "CD68", "Sirpa", "MHCII", "PDL1", "CXCL9")
selected_cts <- c("Alvelor Macrophages", "Interstitial Macrophages", "Monocytes",
                  "Dendritic cells other", "Dendritic cells CD103")

# Run batch and collect stats
stats_table_all <- batch_plot_markers(
  celldata,
  select_ct = selected_cts,
  marker_list = marker_list,
  dir = OUTPUT_PATH,
  image_format = "tiff"
)

# Save all stats
write.csv(stats_table_all, file = file.path(OUTPUT_PATH, "all_markers_stats.csv"), row.names = FALSE)

# Print significant results
sig_stats <- stats_table_all %>% filter(p_value < 0.05)
if (nrow(sig_stats) > 0) {
  cat("Significant differences detected (p < 0.05):\n")
  for (i in seq_len(nrow(sig_stats))) {
    row <- sig_stats[i, ]
    cat(
      paste0("Marker '", row$marker, "' shows significant difference in ",
             row$annotation,
             " (p = ", signif(row$p_value, 3),
             ", ", row$signif, ")\n")
    )
  }
} else {
  cat("No significant differences found (p < 0.05).\n")
}



