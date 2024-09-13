### DE analysis ###

# PACKAGES INSTALLATION
install.packages("readr")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scran")
BiocManager::install("DESeq2")
install.packages("gplots")
install.packages("batchelor")

library(R.utils)
library(readr)
library(BiocManager)
library(SingleCellExperiment)
library(scran)
library(batchelor)
library(gplots)
library(DESeq2)

### LOADING DATA ###

#DIRECTORIES: replace with the appropriate directories

#general files dir
dir_inputfiles <- "./input_files_scrnaseq" #TO ADAPT

#dir_data => where there is the annot_humanAll.csv, features.tsv.gz, matrix.mtx.gz, collectTRI_network.tsv, BP_All_genes.csv, uniprot_ID.tsv, uniprot_to_gene.tab
#dir_output => where we store the outputs of the workflow so we can use to reproduce the results
dir_data <- "./data" #TO ADAPT
dir_output <- "./outputs" #TO ADAPT

#load the pre-processed, down-sampled and normalised SCE object
sce <- readRDS(paste0(directory_output,"/SCE_preproc_norm.rds"))

### IDENTIFICATION OF HVGs, optional step ###

dec <- modelGeneVar(sce,block=sce$cell_type) #block: can be changed
dec <- dec[order(dec$bio, decreasing=TRUE),]

dec2 <- dec[which(dec$FDR <= 0.05 & dec$bio > 0),]

number_hvg <- 2500 #maximal number of HVGs wanted for the downstream analysis, TO ADAPT
prop <- number_hvg/(dim(dec2)[1]) #proportion of genes to report as HVG
hvg <- getTopHVGs(dec2,prop=prop) #genes we will keep as HVG

variablegene <- rep(FALSE,dim(sce)[1])
indexes_hvg <- which(rownames(sce) %in% hvg)
variablegene[indexes_hvg] <- TRUE

jpeg(file=paset0(directory_output,"/variancemean-expression.jpeg"))
plot(dec$mean, dec$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
o <- order(dec$mean)
lines(dec$mean[o], dec$tech[o], col="dodgerblue", lwd=2)
points(dec$mean[variablegene], dec$total[variablegene], type="p", col="red")
dev.off()

sce <- sce[hvg,] #if we want to keep only the HVGs

### BATCH CORRECTION ###
#Adapt the batch(es) to each dataset
sce <- correctExperiments(sce, batch=sce$patient) #TO ADAPT

### SIGNATURE EXTRACTION ###
colLabels(sce) <- sce@colData[[cell_type]]

## FindMarkers Wilcoxon ##

markers_wilcox <- findMarkers(sce,test.type="wilcox",assay.type="logcounts",pval.type="some",min.prop=0.5)

saveRDS(markers_wilcox,paste0(directory_output,"/markers_wilcox.csv"))
markers_wilcox <- readRDS(paste0(directory_output,"/markers_wilcox.csv"))

signature_wilcox <- markers_wilcox
columns_to_keep <- c(1,2,3) #columns from markers_wilcox we want to keep in the signature object
for(i in 1:length(markers_wilcox)){
  signature_wilcox[[i]] <- signature_wilcox[[i]][,columns_to_keep]
}

# #Save the cell type signature
saveRDS(signature_wilcox,paste0(directory_output,"/signature_wilcox.csv"))

#Load the cell type signature
signature_wilcox <- readRDS(paste0(directory_output,"/signature_wilcox.csv"))

## FindMarkers t-test ##

markers_t <- findMarkers(sce,test.type="t",assay.type="logcounts",pval.type="some",min.prop=0.5)

saveRDS(markers_t,paste0(directory_output,"/markers_t.csv"))
markers_t <- readRDS(paste0(directory_output,"/markers_t.csv"))

signature_t <- markers_t
columns_to_keep <- c(1,2,3) #columns from markers_t we want to keep in the signature object
for(i in 1:length(markers_t)){
  signature_t[[i]] <- signature_t[[i]][,columns_to_keep]
}

# #Save the cell type signature
saveRDS(signature_t,paste0(directory_output,"/signature_t.csv"))

#Load the cell type signature
signature_t <- readRDS(paste0(directory_output,"/signature_t.csv"))

## DESeq2 markers ##

colnames(sce@assays@data@listData$counts) <- sce[["cell"]] #TO ADAPT, column name where the cell IDs are stored

#customised sce object to be able to compare 1 VS all cell types (and not 1 VS 1 only)
celltypes <- names(table(sce[[cell_type]]))
res_deseq <- list() #results of the DESeq function
length(res_deseq) <- length(celltypes)
names(res_deseq) <- celltypes

for(i in 1:length(celltypes)){
  type <- celltypes[i]
  custom_sce <- sce
  indexes_other_celltypes <- which(custom_sce[[cell_type]] != type)
  custom_sce[[cell_type]][indexes_other_celltypes] <- "other"
  #Adapt the design argument to each dataset
  sceseq <- DESeqDataSetFromMatrix(countData = custom_sce@assays@data@listData$counts, #TO ADAPT: adapt the design argument to the variables 
                                   colData = custom_sce@colData,
                                   design = ~ patient + annot,
                                   tidy=FALSE)
  #adding of a pseudo count value of 1 to the data 
  sceseq@assays@data@listData$counts <- as.matrix(sceseq@assays@data@listData$counts)+1
  storage.mode(sceseq@assays@data@listData$counts) <- "integer"
  #Adapt the factors to each dataset according to the design argument of the DESeqDataSetFromMatrix function
  sceseq$annot <- factor(sceseq$annot) #TO ADAPT: the variables put in the design argument of the DESeqDataSetFromMatrix function
  sceseq$patient <- factor(sceseq$patient)
  #Adapt the reduced argument to each dataset
  sceseq <- DESeq(sceseq,test="LRT",reduced=~patient,parallel=FALSE,useT=TRUE,minmu=1e-6,minReplicatesForReplace=Inf) #TO ADAPT: adapt the reduced argument to the variables
  res_deseq[i] <- results(sceseq)
}

#Save the DESeq results
saveRDS(res_deseq,paste0(directory_output,"/res_deseq.csv"))

#Load the DESeq results
res_deseq <- readRDS(paste0(directory_output,"/res_deseq.csv"))

#Load the markers_wilcox results
signature_wilcox <- readRDS(paste0(directory_output,"/signature_wilcox.csv"))
signature_deseq <- signature_wilcox #signature_wilcox because it has the right structure, the values will be from res_deseq
columns_to_keep <- c(5,6,2,1,3,4) #columns from res_deseq we want to keep in the signature object, in the right order
for (i in 1:length(signature_wilcox)){ #for each cell type
  res_deseq <- readRDS(paste0(directory_output,"/res_deseq_",i,".csv"))
  signature_deseq[[i]] <- res_deseq[[i]][,columns_to_keep]
}

#Save the DESeq signature
saveRDS(signature_deseq,paste0(directory_output,"/signature_deseq.csv"))

#Load the DESeq signature
signature_deseq <- readRDS(paste0(directory_output,"/signature_deseq.csv"))

## Creation of a common object ##

#put the rownames in the same order in every object (template : first cell type of the wilcoxon test)
signature_all <- SimpleList(signature_wilcox,signature_t,signature_deseq) 
names(signature_all) <- c("wilcoxon","t-test","deseq") #"deseq"
template <- rownames(signature_wilcox[[1]])
for (i in 1:length(signature_all)){ #for each test performed
  for (j in 1:length(signature_all[[i]])){ #for each cell type
    indexes_right_order <- match(template, rownames(signature_all[[i]][[j]]))
    signature_all[[i]][[j]] <- signature_all[[i]][[j]][indexes_right_order,]
  }
}

#Save the ordered signatures
saveRDS(signature_all,paste0(directory_output,"/signature_all.csv"))

#Load the cell signature
signature_all <- readRDS(paste0(directory_output,"/signature_all.csv"))

### SAVING/LOADING THE OBJECT(S) ###
#save the batch-corrected SCE object
saveRDS(sce, paste0(directory_output,"/SCE_batch_corrected.rds"))
#save the HVGs
write.csv(hvg,paste0(directory_output,"/hvg_DEanalysis.csv"))
#save the gene expression signature
saveRDS(signature_all,paste0(directory_output,"/signature_all.csv"))