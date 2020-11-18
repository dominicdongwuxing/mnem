setwd("/Volumes/Macintosh\ HD\ -\ Data/study/yr2sem2/rotation")
library("Matrix")
library("Linnorm")
source('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/code/functions.R')

# raw expression profile data, row: E gene counts; column: individual cell
raw.data <- readMM("GSM2396856_dc_3hr.mtx.txt")

# row: S gene name; corresponding column: names of the cells that were perturbed by this gene
tag <- read.csv(file = "GSM2396856_dc_3hr_cbc_gbc_dict_strict.csv", header = F)
# for faster search, sort all cell names of each sgRNA into a long list
#sorted.cell.names.list <- lapply(tag[,2], function(x) sort(unlist(strsplit(x,','))))
# names of the 32777 cells as columns of raw.data
cell.names <- read.csv(file = "GSM2396856_dc_3hr_cellnames.csv")
# find which sgRNA was detected in raw data by matching cell names to corresponding sgRNA
sgRNA.indices <- lapply(cell.names[,2], function(x) grep(x,tag[,2]))
# convert the S gene indices of cells with no perturbed gene (previously empty) into "$"
sgRNA.indices[unlist(lapply(sgRNA.indices,function(x) length(unlist(x))==0))] <- '$'

# look at some of the properties of data
show.sgRNA.percentage(sgRNA.indices)
show.sgRNA.distribution(sgRNA.indices)

# boolean vector of which valid cells to include (cells with exactly 1 sgRNA match)
sgRNA.inclusion.bool <- unlist(lapply(sgRNA.indices,function(x) length(x) == 1 && x != '$'))
# get which sgRNA is detected in the valid cells 
sgRNA.included <- sgRNA.indices[unlist(sgRNA.inclusion.bool)]

# now convert sgRNA number from sgRNA.included into s gene names 
# by matching this number, which is row number in tag, to a unique S gene
sgene.included <- unlist(lapply(sgRNA.included, function(x) unlist(strsplit(tag[x,1],'_'))[2]))

# get a vector of all 24 S genes, exclude MouseNTC, which is control
sgene.unique <- unique(sgene.included)
sgene.unique <- sgene.unique[sgene.unique != 'MouseNTC']
write.table(sgene.unique,file = 'sgene_unique.txt', sep='\n', row.names = F, col.names = F, quote = F)

# store the s-gene of each included cell into files for later use
write.table(sgene.included,file = "sgene_included.txt", sep='\n', row.names = F, col.names = F, quote = F)
# get and store the s-gene of each included cell, but exclude the control cells
# these are the cells that are in the log odd matrix R
sgene.included.noncontrol <- sgene.included[sgene.included != 'MouseNTC']
write.table(sgene.included.noncontrol,file = "sgene_R.txt", sep='\n', row.names = F, col.names = F, quote = F)

# get a boolean vector of whether each row should be included (median >= 1)
egene.inclusion.bool <- apply(raw.data,1,function(r) median(r) >= 1)
# and use the boolean list to only include rows of raw.data with "True" in boolean list at corresponding position
# also only include columns (cells) with exactly 1 sgRNA matched
included.data <- raw.data[egene.inclusion.bool,sgRNA.inclusion.bool]
# names of the observed 17775 E genes as rows of raw.data
E.genes <- read.csv(file = "GSM2396856_dc_3hr_genenames.csv")
# include useful E genes the same way as did with raw data
egene.included <- E.genes[egene.inclusion.bool,][,2]
# store the valid E-genes for later use
write.table(egene.included,file = "egene_included.txt", sep='\n', row.names = F, col.names = F, quote = F)

# normalize and transform data with Linnorm
normalized.transformed.data <- Linnorm(included.data)
write.table(normalized.transformed.data,'normalized_transofrmed_data.txt',quote = F, sep = '\t', row.names = F, col.names = F)

