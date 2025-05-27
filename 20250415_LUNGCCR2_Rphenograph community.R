
# Rphenograph clustering of normalised, concatenated and scaled data 

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data

## Megan Cole 

#Load packages
library(ggplot2)
library(plyr)
#library(cytofkit)
library(Rphenograph)
print("Phenograph package loaded")

##############################
#### SET GLOBAL VARIABLES ####
##############################
## check the inpt data 
celldistance_calculation<-read.csv("/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile//Results/Spatial_output/20250411_allROIs_25px_HNSCC_neighbours_distanceCalculation_allCells.csv")
if (length(unique(celldistance_calculation$ROI_ID)) == 10) {
  print("可以计算frequence")
} else {
  print("ROI数量不足10个，建议检查数据")
}
dim(celldistance_calculation)

# Name of data file to read in 
path = "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile"

filename = "frequencies_25px_20250411.csv"

experiment_path = "/Results/Spatial_output/"

output_path = "/Results/Spatial_output/"
Sys.setenv(OMP_NUM_THREADS="8") # inside RStudio
Sys.getenv("OMP_NUM_THREADS")
# Load data 
# celldata_fre = read.csv(paste(path, experiment_path, filename, sep = ""))
celldata_fre = read.csv(paste(path, experiment_path, filename, sep = ""), check.names=FALSE)
# celldata_fre = tail(celldata_fre, n = 2000) 
# celldata_fre <- celldata_fre_o[1:2000,]
####################################################################

# Obtain information about the loaded dataset
dim(celldata_fre)
names(celldata_fre)




# Indicate expected length of selectcolumns and select columns for clustering
#Alt MAC	B cell	Cancer	Cl MAC	Cl Mo	DCs cell	Endothelial cell	Int Mo	Mast cell	NK cell	Neutrophils	Non-Cl Mo	T other	Tc	Th	Treg
#check the columns carefully 
relevant_columns = c()
# 13 because note below`
columns_of_interest = 12
selectcolumns = c(grep("B cells", names(celldata_fre), ignore.case = TRUE),
                  grep("Dendritic cells", names(celldata_fre), ignore.case = TRUE), 
                  grep("Endothelium", names(celldata_fre), ignore.case = TRUE),
                  grep("Epithelium", names(celldata_fre), ignore.case = TRUE),
                  grep("Fibroblasts", names(celldata_fre), ignore.case = TRUE),
                  grep("Macrophages", names(celldata_fre), ignore.case = TRUE),
                  grep("Neutrophils", names(celldata_fre), ignore.case = TRUE),
                  grep("T cell CD4", names(celldata_fre), ignore.case = TRUE),
                  grep("T cell CD8", names(celldata_fre), ignore.case = TRUE),
                  grep("Tregs", names(celldata_fre), ignore.case = TRUE), 
                  grep("Tumour", names(celldata_fre), ignore.case = TRUE))
              
# NOTA BENE! WE take both DC classes together with this but this could cause problems down the line. Make sure you select all the right columns here.
                  # grep("Dendritic cells CD103", names(celldata_fre), ignore.case = TRUE))

selectcolumns
length(selectcolumns)
# Check if correct columns have been selected
if(length(selectcolumns) < columns_of_interest){
  print("WARNING: there is a problem with the data selection - missing column name(s)")
} else if(length(selectcolumns) > columns_of_interest){
  # Remove repeat selection of columns 
  selectcolumns = unique(selectcolumns)
  length(unique(selectcolumns))
  
  if(length(selectcolumns) > columns_of_interest){
    # If the length is still above 17, print a warning message that data selection is wrong
    stop(paste("There is a problem wtith data selection - after an attempt to fix, there are still too many columnns selected, expected:", 
               columns_of_interest, "found:", length(selectcolumns)))
  }
  else if(length(selectcolumns) < columns_of_interest){
    stop(paste("There is a problmen with the data selection - missing column name(s), expected:", columns_of_interest,
               "found:", length(selectcolumns)))
  } else{
    print(paste("The list of selected columns has been corrected - there are now", columns_of_interest, "selected columns"))
  }
} else {
  print(paste("The correct number of columns were found - ", columns_of_interest))
}

# Retrieve names of select columns 
names(celldata_fre[,selectcolumns])

# Run Rphenograph on data of select columns
Rphenograph_out <- Rphenograph(celldata_fre[,selectcolumns], k=300)

pheno_temp <- length(as.numeric(membership(Rphenograph_out[[2]])))
pheno_temp
dim(celldata_fre)

celldata_fre$cluster <- as.numeric(membership(Rphenograph_out[[2]]))

# Save data with clustering information added
write.csv(celldata_fre, paste(path, output_path, date,"Rphenograph_HNSCC_output_", max(celldata_fre$cluster), "clusters_k300_", columns_of_interest, "ct_fractions.csv", sep = ""))
print("file saved")


############################################################################################################
