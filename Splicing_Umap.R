

# 加载必要的包
library(jpeg)
library(grid)
# 设置图片文件夹路径
folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/CellCD_markers"

jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)


n_rows <- 4 
n_cols <- 4 

grid.newpage()

pushViewport(viewport(layout = grid.layout(n_rows, n_cols)))


for (i in seq_along(jpeg_files)) {

  img <- readJPEG(jpeg_files[i])
  
  row <- ceiling(i / n_cols)
  col <- ifelse(i %% n_cols == 0, n_cols, i %% n_cols)
  
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  grid.raster(img)
  popViewport()
}

dev.copy(jpeg, "combined_image.jpeg", width = 800, height = 600)  
dev.off()