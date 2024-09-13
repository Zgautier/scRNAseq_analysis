### PRE_PROCESSING AND NORMALISATION OF SC-RNA-SEQ DATASETS ###

### PACKAGES INSTALLATION ###
install.packages("readr")
install.packages("Seurat")
install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scuttle")

library(R.utils)
library(readr)
library(Seurat)
library(BiocManager)
library(SingleCellExperiment)
library(scuttle)
library(biomaRt)

### LOADING DATA ###

#DIRECTORIES: replace with the appropriate directories

#general files dir
dir_inputfiles <- "./input_files_scrnaseq" #TO ADAPT

#dir_data => where there is the annot_humanAll.csv, features.tsv.gz, matrix.mtx.gz, collectTRI_network.tsv, BP_All_genes.csv, uniprot_ID.tsv, uniprot_to_gene.tab
#dir_output => where we store the outputs of the workflow so we can use to reproduce the results
dir_data <- "./data" #TO ADAPT
dir_output <- "./outputs" #TO ADAPT

#with the 3 files matrix.mtx, features.tsv, barcodes.tsv and an annotation matrix
seurdata <- Read10X(dir_data,gene.column=1)
seurobj <- CreateSeuratObject(seurdata)
annot_mat <- read.csv(paste0(dir_data,"/annot_mat.csv")) #load the annotation matrix
colIDcell <- 2 #column of the annotation matrix containing the cell IDs, TO ADAPT
#with a Rdata file (expression matrix)
seurobj <- readRDS(paste0(dir_data,"/Rdata.rds"))
cell_type <- "annot" #name of the annotation column in the object: "cell_type", or "annot"

sce <- as.SingleCellExperiment(seurobj)

#function to replace the " " and "/" by "_" in the cell type names, to do only once and only if necessary
replace_character <- function(sceobj,cell_type){ #cell_type = name of the annotation column in the object
  levels(sceobj[[cell_type]]) <- gsub(" ","_",levels(sceobj[[cell_type]]))
  levels(sceobj[[cell_type]]) <- gsub("/","_",levels(sceobj[[cell_type]]))
  sceobj[[cell_type]] <- gsub(" ","_",sceobj[[cell_type]])
  sceobj[[cell_type]] <- gsub("/","_",sceobj[[cell_type]])
  saveRDS(sceobj,paste0(dir_output,"/renamed_sce.rds"))
  return(sceobj)
}
sce <- replace_character(sce,cell_type) #just done once
#If the genes are in the Ensembl format
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl_genes <- rownames(sce)
gene_ID <- rep(NA, length(ensembl_genes))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                values=ensembl_genes,mart= mart)
genes_to_remove <- which(G_list[,2]=="") #genes for which there is no hgnc symbol
G_list <- G_list[-genes_to_remove,]
genes_in_right_order <- match(rownames(sce),G_list[,1]) #match the indexes of ensembl annotations from G_list to the rownames of sce
genes_in_right_order <- genes_in_right_order[!is.na(genes_in_right_order)] #indexes of G_list corresponding to the rownames of sce
G_list <- G_list[genes_in_right_order,]
indexes_to_keep_sce <- which(rownames(sce) %in% G_list[,1])
sce <- sce[indexes_to_keep_sce,] #keep the genes for which we have a gene hgnc symbol
#save the indexes to keep
saveRDS(indexes_to_keep_sce,paste0(dir_output,"/indexes_to_keep_sce.csv"))
rownames(sce) <- G_list[,2] #replace the rownames by the hgnc symbol

#Save the correctly named (without " " or "/"), and with the gene hgnc IDs and not ensembl IDs
saveRDS(sce,paste0(dir_output,"/renamed_sce.rds"))

#Load the correctly named sce object
sce <- readRDS(paste0(dir_output,"/renamed_sce.rds"))

### QUALITY CONTROL ###

#Mitochondrial transcripts
is.mito <- grepl("MT-",rownames(sce))
summary(is.mito)

sce <- addPerCellQC(sce,subsets=list(Mt=is.mito))

#Library sizes and number of expressed features
libsize.drop <- isOutlier(sce$sum, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$detected, nmads=3, type="lower", log=TRUE)

summary(libsize.drop)
summary(feature.drop)

# save histograms in jpeg format
jpeg(file=paste0(dir_output,"/libsize_hist.jpeg"))
hist(sce$sum, xlab="Library size",
     ylab="Number of cells", breaks=20, main="", col="grey80")
dev.off()
jpeg(file=paste0(dir_output,"/number_features_hist.jpeg"))
hist(sce$detected, xlab="Number of expressed features",
     ylab="Number of cells", breaks=20, main="", col="grey80")
dev.off()
jpeg(file=paste0(dir_output,"/mito_percent_hist.jpeg"))
hist(sce$subsets_Mt_percent, xlab="Mitochondrial proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
dev.off()

mito.drop <- isOutlier(sce$subsets_Mt_percent, nmads=3, type="higher")
summary(mito.drop)

#cell filtering:
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), Remaining=ncol(sce))

### INTERSECTION OF QCs AND ORDER OF THE ANNOTATION MATRIX ###
## When there is an annotation matrix given by the study ##

#enter the column of the annotation matrix where the cell IDs are stored
#suppression of cells filtered in the sce object in the annotation matrix and put in right order the annotation matrix
indexes_annotmat_fromsce <- match(colnames(sce),annot_mat[[colIDcell]]) #indexes of annot_mat corresponding to the cells in sce (or NA if not present)
indexes_annotmat_fromsce <- indexes_annotmat_fromsce[!is.na(indexes_annotmat_fromsce)]
cell_suppressed <- dim(annot_mat)[1]-length(indexes_annotmat_fromsce)
cat('Number of cells suppressed in the annotation matrix because not in sce object: ',
    cell_suppressed, ' cell(s)')

#suppression of cells filtered in the annotation matrix in the sce object
indexes_sce_fromannomat <- match(annot_mat[[colIDcell]],colnames(sce)) #indexes of sce corresponding to the cells in annot_mat (or NA if not present)
indexes_sce_fromannomat <- indexes_sce_fromannomat[!is.na(indexes_sce_fromannomat)]
cell_suppressed <- dim(sce)[2]-length(indexes_sce_fromannomat)
cat('Number of cells suppressed in the sce object because not in the annotation matrix: ',
    cell_suppressed, ' cell(s)')

#Save the index files
write.csv(indexes_annotmat_fromsce,paste0(dir_output,"/indexes_annotmat_fromsce.csv"))
write.csv(indexes_sce_fromannomat,paste0(dir_output,"/indexes_sce_fromannomat.csv"))

#Load the index files
indexes_annotmat_fromsce <- read.csv(paste0(dir_output,"/indexes_annotmat_fromsce.csv"))[-1][[1]]
indexes_sce_fromannomat <- read.csv(paste0(dir_output,"/indexes_sce_fromannomat.csv"))[-1][[1]]

annot_mat <- annot_mat[indexes_annotmat_fromsce,]
sce <- sce[,indexes_sce_fromannomat]

#Following code: needs to be adapted to the available columns of the annotation matrix
#To perform only if there is an annotation matrix in input
coldata <- DataFrame( 
  cell_type=annot_mat[,"cell_type"],  #add all relevant metadata columns
  diagnosis=annot_mat[,"diagnosis"]
)
colData(sce) <- coldata

#Selection of the healthy samples only if necessary
sce <- sce[,which(sce$diagnosis=="Control")]

#Selection of cell types having more than 100 cells
cell_types <- table(sce[[cell_type]])
for (i in 1:length(cell_types)){ #for each cell type
  if (cell_types[[i]] < 100){ #if it has less than 100 cells
    sce <- sce[,which(sce[[cell_type]] != names(cell_types)[i])] #suppress the cells annotated as this cell type
  }
}

#Filtering out low-abundance genes

#Find the smallest cell type, take half its population (N/2) and keep genes expressed in at least N/2 cells in the whole dataset
smaller_pop <- min(table(sce[[cell_type]]))
rownames(sce@assays@data@listData[["counts"]]) <- rownames(sce)
cells_expressing_gene <- rowSums(sce@assays@data@listData[["counts"]] != 0) #for each gene, number of cells in which it is expressed
genes_to_keep <- names(which(cells_expressing_gene>=smaller_pop/2))
cat("number of genes kept after QC:",length(genes_to_keep))

sce <- sce[genes_to_keep,]

### DOWNSAMPLING ON CELL TYPES ###

#Count of each cellular type
celtype <- matrix(ncol=2, nrow= dim(table(sce[[cell_type]])))
colnames(celtype)=c('cell_type','number_cells')
celtype[,1] <- names(table(sce[[cell_type]]))
for(i in 1:dim(table(sce[[cell_type]]))){
  celtype[i,2] <- table(sce[[cell_type]])[[i]]
}
# celltypes_to_remove <- which(as.numeric(celtype[,2])==0)
# celtype <- celtype[-celltypes_to_remove,] #remove cell types that now have 0 counts because had <100 cells 
sceindexes_celtype <- matrix(data=NA,ncol=dim(celtype)[1],nrow=max(as.numeric(celtype[,2])))
colnames(sceindexes_celtype) <- celtype[,1]
for (i in 1:dim(celtype)[1]){ #for each cell type
  size_celltype <- as.numeric(celtype[i,2])
  sceindexes_celtype[1:size_celltype,i] <- which(sce[[cell_type]]==celtype[i,1][[1]])
}

#Save the cell type matrix and indexes of each cell type
write.csv(sceindexes_celtype,paste0(dir_output,"/sceindexes_celtype.csv"))
write.csv(celtype,paste0(dir_output,"/celtype.csv"))

#Load the cell type matrix and indexes of each cell type
sceindexes_celtype <- read.csv(paste0(dir_output,"/sceindexes_celtype.csv"))[,-1]
celtype <- read.csv(paste0(dir_output,"/celtype.csv"))[,-1]

#Downsampling for cellular types
threshold=10000 #maximum wanted population for each cell type, TO ADAPT
sceindexes_sampled <- matrix(data=NA, nrow=threshold, ncol=dim(celtype)[1])
colnames(sceindexes_sampled) <- celtype[,1]
#List of all sampled cell indexes to keep
sampledcells_to_keep <- c()
for (i in 1:dim(celtype)[1]){
  pop_number <- as.numeric(celtype[i,2])
  pop_indexes <- sceindexes_celtype[1:pop_number,i]
  if (pop_number >= threshold){
    sample_indexes <- sample(pop_indexes,threshold,replace=FALSE)
    sceindexes_sampled[,i] <- sample_indexes
    sampledcells_to_keep <- append(sampledcells_to_keep, sample_indexes)
  }
  else{
    sceindexes_sampled[1:pop_number,i] <- sceindexes_celtype[1:pop_number,i]
    sampledcells_to_keep <- append(sampledcells_to_keep,sceindexes_celtype[1:pop_number,i])
  }
}

#Save the cell indexes to keep in each cell type and the total cell indexes to keep
write.csv(sampledcells_to_keep,paste0(dir_output,"/sampledcells_to_keep.csv"))

#Load the cell indexes to keep in each cell type and the total cell indexes to keep
sampledcells_to_keep <- read.csv(paste0(dir_output,"/sampledcells_to_keep.csv"))[,-1]

cells_to_keep <- colnames(sce)[sampledcells_to_keep]

#Save the cell IDs to keep
write.csv(cells_to_keep,paste0(dir_output,"/cells_to_keep.csv"))

sce <- sce[,sampledcells_to_keep]
cat('Number of cells in the sampled matrix: ', dim(sce)[2], ' cells')

#Save the pre-processed and non-normalised sce object
saveRDS(sce,paste0(dir_output,"/preprocessed_sce.rds"))

#Save the pre-processed and non-normalised sce object
sce <- readRDS(paste0(dir_output,"/preprocessed_sce.rds"))

### NORMALISATION OF CELL-SPECIFIC BIASES ###
sce <- scuttle::logNormCounts((sce), librarySizeFactors((sce)))

### SAVING/LOADING THE OBJECT(S) ###
#save the pre-processed, down-sampled and normalised SCE object
saveRDS(sce, paste0(dir_output,"/SCE_preproc_norm.rds"))
