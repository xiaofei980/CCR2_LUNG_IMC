
### code one : slow  but works   ,and can load many images 

library(jpeg)
library(grid)

folder_path <-  "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/TUMCD_markers"

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

## code two :   work  and fast 

library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/MacDendCD_markers"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)

n_rows <- 3
n_cols <- 4


images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}
#################### Extra code for splicing ,better format 
library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/MacDendCD_markers/For_splicing"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)


n_rows <- 1
n_cols <- 3

images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}


####################CellCD

library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/CellCD_markers"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)

n_rows <- 3
n_cols <- 4


images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}
#################### Extra code for splicing ,better format 
library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/CellCD_markers/For_splicing"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)


n_rows <- 1
n_cols <- 4

images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}

########MacDendTumCD

library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/MacDendTumCD_markers"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)

n_rows <- 3
n_cols <- 3


images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}
#################### Extra code for splicing ,better format 
library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/MacDendTumCD_markers/For_splicing"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)


n_rows <- 3
n_cols <- 1

images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}

##### T cell cd

library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/Results/Figures_output/UMAP/T_Cell_markers"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)

n_rows <- 3
n_cols <- 4


images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}
######## celldatavs map
library(magick)

gc()  # 触发垃圾回收，释放内存
rm(list = ls())  # 清空所有变量

folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/240920_Sahai_CCR2_IMC/Projectfile/map_rawImage/EPTumour"
jpeg_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)

n_rows <- 2
n_cols <- 5


images_per_combined <- n_rows * n_cols  # 每张组合图包含的图片数量
num_batches <- ceiling(length(jpeg_files) / images_per_combined)  # 计算需要生成多少张组合图

for (i in 1:num_batches) {
  start_idx <- (i - 1) * images_per_combined + 1
  end_idx <- min(i * images_per_combined, length(jpeg_files))
  
  batch_files <- jpeg_files[start_idx:end_idx]  # 获取当前批次的图片
  
  images <- lapply(batch_files, image_read)  # 读取当前批次的图片
  
  combined_image <- image_montage(
    image_join(images),
    tile = paste0(n_cols, "x", n_rows),
    geometry = "x500",
    bg = "white"
  )
  
  output_path <- file.path(folder_path, paste0("combined_image_", i, ".jpg"))
  image_write(combined_image, path = output_path, format = "jpg")
  cat("Image saved at:", output_path, "\n")
}


