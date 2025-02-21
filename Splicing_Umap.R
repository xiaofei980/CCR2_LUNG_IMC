library(magick)

# 设置图片所在的文件夹路径
folder_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Results/Annotation_output/UMAP/T_Cell_markers"

# 读取文件夹中的所有 JPEG 图片
image_files <- list.files(folder_path, pattern = "\\.jpeg$", full.names = TRUE)

# 选出包含关键字的图片（"group", "cluster", "annotation"）
priority_images <- grep("group|cluster|annotation", image_files, value = TRUE, ignore.case = TRUE)

# 其余普通图片
remaining_images <- setdiff(image_files, priority_images)

# 组合所有图片，确保关键字图片排在前面
ordered_images <- c(priority_images, remaining_images)

# 设定拼接的行数和列数
nrow <- 4  # 设定行数
ncol <- 4  # 设定列数

# 读取所有图像
image_list <- lapply(ordered_images, image_read)

# 确保图像数量和网格匹配（如果不足，可以填充空白）
num_images <- length(image_list)
grid_size <- nrow * ncol

if (num_images < grid_size) {
  empty_image <- image_blank(width = image_info(image_list[[1]])$width, 
                             height = image_info(image_list[[1]])$height, 
                             color = "white")  # 创建空白占位符
  image_list <- c(image_list, rep(list(empty_image), grid_size - num_images))
}

# 使用 image_montage 进行网格拼接
combined_image <- image_montage(image_join(image_list), tile = paste0(ncol, "x", nrow), geometry = "100%x100%")

# 显示拼接结果
print(combined_image)

# 保存最终大图
output_path <- file.path(folder_path, "combined_image.jpeg")
image_write(combined_image, path = output_path, format = "jpeg")
