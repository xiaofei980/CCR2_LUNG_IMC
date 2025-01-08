# imagemtadata is for getting the width and Height  for the images
library(exiftoolr)
library(stringr)

# Script used to extract filename, height and width of images
base = "/mnt/Data1/groupfebe/runs/Xiaofei/CCR2_LUNG/Data/"
image_path = file.path(base,"Normalization_input/full_stack.ome")
OUTPUT_path = file.path(base, "Normalization_input/imagemetadata.csv")


image_metadata = exif_read(path = file.path(image_path, list.files(image_path)))
image_metadata = image_metadata[,c("FileName", "ImageWidth", "ImageHeight")]
colnames(image_metadata) = c("filename", "width", "height")
image_metadata$filename = str_remove(string = image_metadata$filename,
                                     pattern = "_full_stack.ome.tiff")
View(image_metadata)
write.csv(image_metadata, file = OUTPUT_path, row.names = F)




