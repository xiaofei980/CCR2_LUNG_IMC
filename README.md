# CCR2_LUNG_IMC-
The code is for  IMC data from CCR2KO lung samples
step1 : get reorganize data from segmentation folder : fullstack.ome--> imagematadata(include the width and length for later using )  and Allcellmask (include all the mean intensity and domain information)
step2:  runs normalization(in this step mainly consider those markers which performs not good ,  for example LAG3 ,TIM3 etc,if not expressed ,the mean intensity can be upregulating.)  and tSNE (depend on the data sized )and Rphenotype(consider using unsupervised or superviese cluster method ,which is better , in this script ,i use unsupervised RPhenograph) 
Annotation : plotscatter to set up the threshold , umpa visulaization checking , cross checking with FIJI . 
