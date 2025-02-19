## Sampling 100,000 data for generating tSNE, this is important if the all dataset is large , running the whole data tSNE will cause the full C stack . 
#For my experience, 200k is ok but 1m will definetely collapse.
print(Sys.time())
print(paste("C stack size:", Cstack_info()["size"]))
print("Loading packages...")
library(stringr)
library(dplyr)
library(Rtsne)
print("Packages loaded")


# Path for loading data and saving results 
INPUT_PATH = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Data/tSNE_input/"
OUTPUT_PATH = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/tSNE_output/"

print(Sys.time())
print("Loading celldata...")
# Load celldata excluding tSNE coordinates
celldata = read.csv(file.path(INPUT_PATH, "250218_normalised_concatenated_scaled_data_001threshold.csv"),
                    stringsAsFactors = F)


print("Celldata loaded")
## in this dataset ,i will not set sample size due to in total 16w ,similar to 10w 
# # Sample size for t-SNE analysis
# sample_size <- 100000  # Adjust this according t your dataset size and computational resources1
# 
# # Sample data at the beginning
# celldata_sampled <- celldata[sample(nrow(celldata), sample_size), ]
# 
# # Select columns to be used for t-SNE generation:
# # All markers except BG, DNA and domain
# tSNEcolumns = grep("MI_", names(celldata_sampled), value = TRUE)
# exclCOLs = c("MI_Argon80", "MI_Xenon131", "MI_Xenon132", "MI_RuO4", "MI_DNA1",
#              "MI_DNA2", "MI_TumourminusInterface", "MI_NormalminusInterface",
#              "MI_InterfaceImage")
# tSNEcolumns = tSNEcolumns[!tSNEcolumns %in% exclCOLs]

print(Sys.time())

# # All markers except BG, DNA and domain
 tSNEcolumns = grep("MI_", names(celldata), value = TRUE)
 exclCOLs = c("MI_Argon80", "MI_Xenon131", "MI_Xenon132", "MI_RuO4", "MI_DNA1",
            "MI_DNA2", "MI_TumourminusInterface", "MI_NormalminusInterface",
           "MI_InterfaceImage")
 tSNEcolumns = tSNEcolumns[!tSNEcolumns %in% exclCOLs]


## perplexity from 100 to 60 , Max_iter from 1500 to 1200,  due to C stack limit
tSNE_results = Rtsne(celldata[,tSNEcolumns], num_threads = 4,
                     perplexity = 100, verbose = TRUE,
                     check_duplicates = FALSE, max_iter = 1500)

print("Running t-SNE")


print("t-SNE completed")

tSNE1 = tSNE_results$Y[,1]
tSNE2 = tSNE_results$Y[,2]

tSNE = cbind(tSNE1, tSNE2)



print(Sys.time())
print("Finished!")
date <- str_sub(gsub("-", "", Sys.Date()), -6, -1)
#output_file <- paste0(OUTPUT_PATH, date, "_celldata.csv")
output_file2 <- paste0(OUTPUT_PATH, date, "_tSNE.csv")
# Write the data to the specified location
#write.csv(celldata_sampled, file = output_file, row.names = FALSE)
write.csv(tSNE,file=output_file2,row.names = FALSE)
# Print out the full path of the output file
print(paste("File saved at:", output_file2))

print(Sys.time())
print("Finished!")


