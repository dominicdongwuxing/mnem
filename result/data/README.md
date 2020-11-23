# Processed data

This directory contains the processed data of the dendritic cell from the Perturb-seq paper (Dixit et al., 2018). A link to the paper can be found here (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5181115/).

To obtain "normalized_transformed_data.csv", the preprocessData function is used, which converts raw data matrix into preprocessed data matrix by selecting cells with only one sgRNA detected, and deleting E genes (rows) with median < 1, and then perform normalization and transformation of the data, using the Linnorm package. 

 The function uses "GSM2396856_dc_3hr.mtx.txt" (the gene expression count data matrix); "GSM2396856_dc_3hr_cbc_gbc_dict_strict.csv" (lists for each sgRNA, the cell names where it is detected) and "GSM2396856_dc_3hr_cellnames.csv" (the cell names for the data matrix), which can to be downloaded from single cell protal (https://singlecell.broadinstitute.org/single_cell/study/SCP24/perturb-seq#study-download) and need to be stored in the working directory. 
 
 To obtain "loggodd.csv", the function computeLogodd is used, which uses the preprocessed data matrix from above, and convert the normalized transformed data matrix into log odd matrix. 
 