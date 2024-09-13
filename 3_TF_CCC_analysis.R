### CCC ###

# PACKAGES INSTALLATION
install.packages("readr")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
install.packages("gplots")
install.packages("batchelor")
BiocManager::install("viper")
BiocManager::install("enrichR")
BiocManager::install("OmnipathR")

library(R.utils)
library(readr)
library(BiocManager)
library(SingleCellExperiment)
library(batchelor)
library(gplots)
library(viper)
library(stringr)
library(tibble)
library(purrr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(magrittr)
library(liana)
library(ggplot2)
library(patchwork)
library(knitr)
library(enrichR)

### LOADING DATA ###

#DIRECTORIES: replace with the appropriate directories

#general files dir
dir_inputfiles <- "./input_files_scrnaseq" #TO ADAPT

#dir_data => where there is the annot_humanAll.csv, features.tsv.gz, matrix.mtx.gz, collectTRI_network.tsv, BP_All_genes.csv, uniprot_ID.tsv, uniprot_to_gene.tab
#dir_output => where we store the outputs of the workflow so we can use to reproduce the results
dir_data <- "./data" #TO ADAPT
dir_output <- "./outputs" #TO ADAPT

#load the pre-processed, down-sampled, normalised and batch-corrected SCE object
sce <- readRDS(paste0(directory_output,"/SCE_batch_corrected.rds"))
#load the HVGs
hvg <- read_csv(paste0(directory_output,"/hvg_DEanalysis.csv"))[[2]]
#load the gene expression signature
signature_all <- readRDS(paste0(directory_output,"/signature_all.csv"))

### SELECTION OF THE DEGs ###

#Combining the 4 tests for the DEGs

#choose one of the signatures
signature_all[["common"]] <- signature_all$`t-test`

for(i in 1:length(signature_all$`t-test`)){ 
  colnames(signature_all[["common"]][[i]]) <- c("pvalue", "FDR", "LFC")
}

for(celltype in names(signature_all$`t-test`)){ #for each cell type
  pvalue_threshold <- 0.05
  #Wilcoxon
  selected_wilcox <- which(signature_all$wilcoxon[[celltype]]$p.value <= pvalue_threshold)
  #t-test
  selected_t <- which(signature_all$`t-test`[[celltype]]$p.value <= pvalue_threshold)
  #deseq test
  selected_deseq <- which(signature_all$deseq[[celltype]]$pvalue <= pvalue_threshold)
  
  selected_genes <- intersect(selected_t,intersect(selected_wilcox,selected_deseq))
  
  signature_all$common[[celltype]] <- signature_all$common[[celltype]][selected_genes,]
}

signature <- signature_all$common

### TF ACTIVITY ###

#CollectTRI-derived regulons

df <- read.delim(paste0(directory_inputfiles,"/collectTRI_network.tsv"))

df2regulon <- function(df) { #dataframes to regulon objects
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}

regulon_collectri <- df2regulon(df)

#msVIPER

list_mrs <- list() #list of msVIPER results for each cell type
list_mrs_df <- list() #list of msVIPER results for each cell type in a dataframe

for(i in 1:length(signature)){ #for each cell type, analysis with all genes
  signature_celltype <- signature[[i]][,2] #FDR column
  sig_logFoldChanges <- signature[[i]][,3] #logFoldChanges column
  signature_celltype <- qnorm(signature_celltype/2, lower.tail=FALSE)*sign(sig_logFoldChanges)
  names(signature_celltype) <- rownames(signature[[i]])
  mrs <- msviper(ges = signature_celltype, regulon = regulon_collectri, ges.filter = F, minsize = 4,verbose=F)
  mrs_df <- data.frame('size' = mrs$es$size,
                       'p.value' = mrs$es$p.value,
                       'fdr' = p.adjust(mrs$es$p.value, method = 'fdr'),
                       'nes' = mrs$es$nes)
  list_mrs[[i]] <- mrs
  list_mrs_df[[i]] <- mrs_df
  names(list_mrs)[i] <- names(signature)[i]
  names(list_mrs_df)[i] <- names(signature)[i]
}

#Save the msVIPER output
saveRDS(list_mrs,paste0(directory_output,"/list_mrs.csv"))

#Save the msVIPER output in DataFrames
saveRDS(list_mrs_df,paste0(directory_output,"/list_mrs_df.csv"))

#Load the msVIPER output
list_mrs <- readRDS(paste0(directory_output,"/list_mrs.csv"))

#Load the msVIPER output in DataFrames
list_mrs_df <- readRDS(paste0(directory_output,"/list_mrs_df.csv"))

### CELL CELL COMMUNICATION ###

#selection of DEGs in the signature (adaptive threshold depending on the cell type)
dim_deg <- c()
min <- 0.15 #minimum logFC threshold to define a DEG
max <- 2.5 #maximum logFC threshold to define a DEG

for (i in 1:length(signature)){ #for each cell type
  threshold <- min
  
  #with the min threshold
  indexes_logFC <- which(signature[[i]]$LFC>threshold) #indexes of genes having a logFC>threshold
  selected_genes <- signature[[i]][indexes_logFC,] #selection of those genes
  
  while( (threshold <= max) && (dim(selected_genes)[1] >= 500) ){
    #while the logFC threshold is <= max and the length of DEGs >= 500
    threshold <- threshold + 0.01
    indexes_logFC <- which(signature[[i]]$LFC>threshold) #indexes of genes having a logFC>threshold
    selected_genes <- signature[[i]][indexes_logFC,] #selection of those genes
  }
  
  celltype <- names(signature)[i]
  dir_file <- paste0(directory_data,"/Table_DEGs_t/Signalling_DEGs_",celltype,".csv")
  if (!(file.exists(dir_file))){
    dir.create(file.path(directory_data, "Table_DEGs_t"))
  }
  table_deg <- as.data.frame(selected_genes)
  table_deg <- cbind(table_deg, gene = rownames(selected_genes))
  saveRDS(table_deg, dir_file) #save the results as a dataframe in the Table_DEGs file
  dim_deg <- append(dim_deg,dim(table_deg)[1])
}

hist(dim_deg,breaks=20, xlab="number of DEGs per cell type with T-test, 0.15 min threhsold")

colLabels(sce) <- sce@colData[[cell_type]] #cell type annotation as labels

#workflow on the LIANA tool, aggregating the results of the different databases used
liana_results <- liana_wrap(sce)
aggregated_liana <- liana_aggregate(liana_results)
aggregated_liana <- cbind(aggregated_liana, tool = rep("LIANA",dim(aggregated_liana)[1]))
LIANA_net <- aggregated_liana

#Save the LIANA results
saveRDS(liana_results, paste0(directory_output,"/liana_results.csv"))
saveRDS(LIANA_net,paste0(directory_output,"/LIANA_net.csv")) #we choose to use the unfiltered object with all interactions

#Load the LIANA results
liana_results <- readRDS(paste0(directory_output,"/liana_results.csv"))
LIANA_net <- readRDS(paste0(directory_output,"/LIANA_net.csv"))

#0. Merge the LIANA results
LIANA_net <- readRDS(paste0(directory_output,"/LIANA_net.csv"))
LIANA_net <- LIANA_net %>%
  mutate(interaction_pair = receptor.complex) %>% #create new column interaction_pair=receptor complex, nb(rows) not affected
  separate_rows(receptor.complex, sep = "_") #separate in different rows (receptor complex -> if several receptors, separated in several rows, interaction pair are kept in the interaction pair column)
LIANA_net <- LIANA_net[,c(1:4,17,14)] #selection of source, target, ligand_complex, receptor, interaction_pair, sca_rank
colnames(LIANA_net)[c(3,4,6)] <- c("ligand","receptor","score_interaction")
LIANA_CCC_merged_results <- LIANA_net
saveRDS(LIANA_CCC_merged_results, paste0(directory_output,"/LIANA_CCC_merged_results.csv"))

#1. Get the receptors from LIANA DEGs for each cell type

load_tables <- function(celltype, ccc_file){
  #celltype = name of the cell type
  #ccc_file = output of "0_Getting_receptors_from_LIANA.R", named "LIANA_CCC_merged_results.csv"
  ###load tables
  CCC_NH_results <- as.data.frame(ccc_file)
  deg_dir <- paste0(directory_data,"/Table_DEGs/Signalling_DEGs_",celltype,".csv")
  deg_file <- readRDS(deg_dir)
  ### extract data for the given cell type
  CCC_rec <- CCC_NH_results %>% filter(target == celltype) %>% #take the interactions for which the target is the cell type of interest
    dplyr::select(source,target,ligand,receptor,score_interaction) #selection of the following info: source,target,ligand,receptor,score_interaction,tool
  
  #selection of receptors in the DEGs file
  DEG_receptors <- deg_file %>%
    filter(gene %in% CCC_rec$receptor) #with logFC or AUC
  # Save the intersection of receptors and DEGs for each cell type
  file_name <- paste0(directory_data,"/Table_receptors/", celltype, "_DEG_receptors.csv")
  if (!(file.exists(file_name))){
    dir.create(file.path(directory_data, "Table_receptors"))
  }
  saveRDS(DEG_receptors, file_name) #write the result in the Table_receptors file
}
## load the DEGs for every cell type
LIANA_CCC_merged_results <- readRDS(paste0(directory_output,"/LIANA_CCC_merged_results.csv"))
signature_all <- readRDS(paste0(directory_output,"/signature_all.csv"))

for(i in 1:length(signature)){ #for each cell type
  celltype <- names(signature)[i]
  load_tables(celltype, LIANA_CCC_merged_results)
}

#2. Receptor-TF downstream pathway
BP_Genes_GO_BP <- read.delim(paste0(directory_inputfiles,"/BP_All_genes.csv"),sep=",")

Reactome_2024 <- read.delim(paste0(directory_inputfiles,"/reactome_2024.tsv"))
colnames(Reactome_2024) <- NULL
Reactome <- tibble()
colnames(Reactome) <- c("SYMBOL", "R_HSA_Terms")
for(i in 1:length(Reactome_2024[[1]])){
  reactome_line <- Reactome_2024[[1]][i]
  term <- str_split(reactome_line,fixed("\t"))[[1]][1]
  R_HSA_Terms <- word(term,-1) #selection of the R-HSA term only
  SYMBOL <- str_split(reactome_line,fixed("\t"))[[1]][2] #selection of the gene name
  frame <- data.frame(SYMBOL,R_HSA_Terms)
  Reactome <- bind_rows(Reactome,frame)
}
#Save the gene symbol to Reactome R-HSA terms table
write_csv(Reactome, paste0(directory_inputfiles,"/Reactome_RHSA_2024.csv"))

#Load the gene symbol to Reactome R-HSa terms table
Reactome <- read.delim(paste0(directory_inputfiles,"/Reactome_RHSA_2024.csv"),sep=",")

#BP from GO
# Define function to get enriched genes from receptors and transcription factors (TFs)
get_enriched_genes_GO <- function(receptor_file, tf_file, output_file, output_file2) {
  
  #receptor_file = output from the receptor analysis performed from the DEGs
  #tf_file = output from the CollecTRI analysis for the cell type of interest
  #output_file = directory of the output "DEGs_filt_by_path/%s_GO_BP_23_enriched_rec_protein.csv"
  #output_file2 = directory of the output "Enrich_R_tables/%s_GO_BP_23_TF_receptors_rec_protein.csv"
  
  # Load receptor data
  #gene replaced by "receptor" in my case
  receptors <- receptor_file %>%
    dplyr::select(gene, LFC) %>%
    dplyr::rename(SYMBOL = gene, RWRS = LFC)
  
  # Load TF data and process
  #score replaced by nes = normalised enrichment score, in my case
  tfs <- tf_file[which(tf_file$fdr<=0.05),] %>%
    dplyr::select(source, nes) %>%
    dplyr::arrange(desc(nes)) %>% #we take only the positive ones to have the upregulated TFs
    dplyr::slice(1:10) %>% #we take the top 10 up-regulated TFs
    dplyr::rename(SYMBOL = source, RWRS = nes)

  # Combine receptor and TF data
  rec_tf <- dplyr::full_join(receptors, tfs) #concatenate the 2
  
  # Perform enrichment analysis
  seed <- enrichr(rec_tf$SYMBOL, "GO_Biological_Process_2023")
  seed_tab <- as.data.frame(seed)
  
  # Filter for significant results
  significant_genes <- seed_tab %>%
    dplyr::filter(GO_Biological_Process_2023.Adjusted.P.value < 0.05) %>%
    dplyr::pull(GO_Biological_Process_2023.Term) # Extracting just the Term column
  go_terms <- gsub(".*\\(GO:(\\d+)\\).*", "GO:\\1", significant_genes) #selection of the R-HSA-terms only
  
  significant_genes <- BP_Genes_GO_BP %>%
    dplyr::filter(GO_Terms %in% go_terms) #selection of the lines of the BP_Genes_GO_BP database corresponding to the GO-terms selected
  
  saveRDS(significant_genes, output_file)
  saveRDS(seed_tab, output_file2)
}

#TF activity output
list_mrs_df <- readRDS(paste0(directory_output,"/list_mrs_df.csv"))
#Receptor directory
rec_dir <- paste0(directory_data,"/Table_receptors/")

#Process for each cell type
for(i in 1:length(list_mrs_df)){ #for each cell type
  #Receptor file input
  celltype <- names(list_mrs_df)[i]
  rec_dir_for_this_celltype <- paste0(rec_dir, celltype,"_DEG_receptors.csv")
  receptor_file <- readRDS(rec_dir_for_this_celltype)
  
  #TF file input
  tf_file <- as.data.frame(list_mrs_df[[i]])
  tf_file <- cbind(tf_file, source = rownames(tf_file)) #add the rownames (TFs) as a column
  
  #Output file inputs
  output_file <- paste0(directory_data,"/DEGs_filt_by_path/",celltype,"_GO_BP_23_enriched_rec_protein.csv")
  output_file2 <- paste0(directory_data,"/Enrich_R_tables/",celltype,"_GO_BP_23_TF_receptors_rec_protein.csv")
  
  if (!(file.exists(output_file))){
    dir.create(file.path(directory_data,"DEGs_filt_by_path"))
  }
  if (!(file.exists(output_file2))){
    dir.create(file.path(directory_data,"Enrich_R_tables"))
  }
  
  get_enriched_genes_GO(receptor_file, tf_file, output_file, output_file2)
}

#Reactome DB
# Define function to get enriched genes from receptors and transcription factors (TFs)
get_enriched_genes_Reactome <- function(receptor_file,tf_file,reactome_outputfile) {
  
  #receptor_file = output from the receptor analysis performed from the DEGs
  #tf_file = output from the CollecTRI analysis for the cell type of interest
  #output_file = directory of the output "DEGs_filt_by_path/%s_GO_BP_23_enriched_rec_protein.csv"
  #output_file2 = directory of the output "Enrich_R_tables/%s_GO_BP_23_TF_receptors_rec_protein.csv"
  
  # Load receptor data
  #gene replaced by "receptor" in my case
  receptors <- receptor_file %>%
    dplyr::select(gene, LFC) %>%
    rename(SYMBOL = gene, RWRS = LFC)
  
  # Load TF data and process
  #score replaced by nes = normalised enrichment score, in my case
  tfs <- tf_file[which(tf_file$fdr<=0.05),] %>%
    dplyr::select(source, nes) %>%
    arrange(desc(nes)) %>% #we take only the positive ones to have the upregulated TFs
    dplyr::slice(1:10) %>% #we take the top 10 up-regulated TFs
    rename(SYMBOL = source, RWRS = nes)

  # Combine receptor and TF data
  rec_tf <- full_join(receptors, tfs) #concatenate the 2
  
  # Perform enrichment analysis
  seed <- enrichr(rec_tf$SYMBOL, "Reactome_2022")
  significant_genes <- as.data.frame(seed) %>%
    filter(Reactome_2022.Adjusted.P.value < 0.05) %>%
    pull(Reactome_2022.Term) # Extracting just the Term column
  rhsa_terms <- word(significant_genes, -1) #selection of the RHSA-terms only
  write_csv(as.data.frame(rhsa_terms), reactome_outputfile)
}

specific_terms <- function(rhsa_specific_terms,output_file,output_file_2){ #with the selection of the most specific Reactome terms only
  
  significant_genes <- Reactome %>%
    filter(R_HSA_term %in% rhsa_specific_terms[[1]])  #selection of the lines of the Reactome database corresponding to the terms selected
  
  print("significant genes :")
  print(length(unique(significant_genes$SYMBOL)))
  write_csv(significant_genes, output_file)
  write_csv(seed_tab, output_file_2)
}

not_specific_terms <- function(rhsa_not_specific_terms,output_file,output_file_2){ #All reactome terms selected
  
  significant_genes <- Reactome %>%
    filter(R_HSA_term %in% rhsa_not_specific_terms[[1]])  #selection of the lines of the Reactome database corresponding to the GO-terms selected
  
  print("significant genes :")
  print(length(unique(significant_genes$SYMBOL)))
  write_csv(significant_genes, output_file)
  write_csv(seed_tab, output_file_2)
}

#TF activity output
list_mrs_df <- readRDS(paste0(directory_output,"/list_mrs_df.csv"))
#Receptor directory
rec_dir <- paste0(directory_data,"/Table_receptors/")
#Specific Reactome terms directory
reactome_dir_specific <- paste0(directory_output,"/filtered_reactome_terms/RHSA_") #specific terms
reactome_dir_not_specific <- paste0(directory_output,"/reactome_terms/RHSA_") #not specific terms

#Process for each cell type
for(i in 1:length(list_mrs_df)){ #for each cell type
  #Receptor file input
  celltype <- names(list_mrs_df)[i]
  rec_dir_for_this_celltype <- paste0(rec_dir, celltype,"_DEG_receptors.csv")
  receptor_file <- readRDS(rec_dir_for_this_celltype)
  
  #TF file input
  tf_file <- as.data.frame(list_mrs_df[[i]])
  tf_file <- cbind(tf_file, source = rownames(tf_file)) #add the rownames (TFs) as a column
  
  #Output file inputs
  output_file <- paste0(directory_data,"/DEGs_filt_by_path/",celltype,"_Reactome_24_enriched_rec_protein.csv")
  output_file_2 <- paste0(directory_data,"/Enrich_R_tables/",celltype,"_Reactome_24_receptors_rec_protein.csv")
  reactome_outputfile <- paste0(directory_output,"/reactome_terms/RHSA_",celltype)
  
  if (!(file.exists(output_file))){
    dir.create(file.path(directory_data,"DEGs_filt_by_path"))
  }
  if (!(file.exists(output_file_2))){
    dir.create(file.path(directory_data,"Enrich_R_tables"))
  }

  get_enriched_genes_Reactome(receptor_file, tf_file,reactome_outputfile)
  
  # #Specific Reactome terms 
  # rhsa_specific_terms <- read_csv(paste0(reactome_dir_specific,celltype))
  # specific_terms(rhsa_specific_terms,output_file,output_file_2)
  
  #Not Specific Reactome terms
  rhsa_not_specific_terms <- read_csv(paste0(reactome_dir_not_specific,celltype))
  not_specific_terms(rhsa_not_specific_terms,output_file,output_file_2)
}

#3. Get the seed nodes input for sc-phuEGO

process_cell_type <- function(celltype){
  
  # Directories of the input files
  degs_file <- paste0(directory_data,"/Table_DEGs/Signalling_DEGs_",celltype,".csv")
  # significantgenes_file <- paste0(directory_data,"/DEGs_filt_by_path/",celltype,"_Reactome_24_spec_enriched_rec_protein.csv")
  significantgenes_file <- paste0(directory_data,"/DEGs_filt_by_path/",celltype,"_Reactome_24_enriched_rec_protein.csv")
  # significantgenes_file <- paste0(directory_data,"/DEGs_filt_by_path/",celltype,"_GO_BP_23_enriched_rec_protein.csv")
  tfs_file <- paste0(directory_output,"/list_mrs_df.csv")
  receptor_file <- paste0(directory_data,"/Table_receptors/",celltype,"_DEG_receptors.csv")
  
  number_genetypes <- List()
  number_degs <- c()
  number_derec <- c()
  number_tfs <- c()
  
  #Load data
  uniprot_to_gene <- read.table(paste0(directory_inputfiles,"/uniprot_to_gene.tab"), header = F, sep = "\t") %>%
    dplyr::rename(Protein = V1, SYMBOL = V2) %>%
    #rename(V1 = "Protein", V2 = "SYMBOL") %>% 
    dplyr::distinct()  # Ensure no duplicate entries
  UNIPROT_ID <- read.delim(file = paste0(directory_inputfiles,"/uniprot_ID.tsv")) %>%
    dplyr::rename(Protein = Entry) %>%
    #rename(Entry = "Protein") %>%
    dplyr::distinct()  # Ensure no duplicate entries
  uniprot_to_gene <- dplyr::right_join(uniprot_to_gene, UNIPROT_ID, by = "Protein") #table with Uniprot IDs and gene names
  
  # Ensure uniprot_to_gene has unique SYMBOL entries
  uniprot_to_gene <- uniprot_to_gene %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(Protein = dplyr::first(Protein))  # Resolve many-to-many relationships
  
  # Load DEGs
  degs_data <- readRDS(degs_file) %>%
    dplyr::select(SYMBOL=gene,RWRS=LFC)  #DEGs for each cell type
  # Load significant genes
  sig_genes_data <- read.delim(significantgenes_file,sep=",") %>%
    dplyr::select(SYMBOL = SYMBOL)  #DEGs for each cell type
  # sig_genes_data <- readRDS(significantgenes_file) %>%
  #   dplyr::select(SYMBOL = SYMBOL)  #DEGs for each cell type
  
  # Take the intersection
  degs_filtered <- unique(dplyr::inner_join(degs_data, sig_genes_data, by = "SYMBOL"))
  number_degs <- dim(degs_filtered)[1]
  
  # Add receptors + TFs to create the seeds
  #TF file input
  list_mrs_df <- readRDS(tfs_file)
  tfs <- as.data.frame(list_mrs_df[[celltype]])
  tfs <- cbind(tfs,SYMBOL = rownames(tfs)) %>% #add the rownames (TF gene name) as a column
    dplyr::arrange(desc(nes)) %>% #we take only the positive ones to have the upregulated TFs
    dplyr::slice(1:10) %>% #we take the top 10 up-regulated TFs
    dplyr::select(SYMBOL,nes) %>%
    dplyr::rename(RWRS = nes)
  
  number_tfs <- dim(tfs)[1]
  
  #receptor file input
  receptors <- readRDS(receptor_file) %>%
    dplyr::select(gene, LFC) %>%
    dplyr::rename(SYMBOL = gene, RWRS = LFC)
  number_derec <- dim(receptors)[1]
  
  number_genetypes[[celltype]] <- c(number_degs,number_tfs,number_derec)
  
  #join receptors and TFs 
  rec_tf <- rbind(receptors, tfs)
  
  # Merge tables
  merged_data <- dplyr::full_join(degs_filtered, rec_tf,relationship = "many-to-many") %>%
    dplyr::inner_join(uniprot_to_gene, by = "SYMBOL") %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  output_file <- paste0(directory_data,"/Table_significant_DEGs/",celltype,"_significant_DEGs.csv")
  # Save the data
  output_file_csv <- paste0(directory_data,"/seed_nodes_scphuEGOCSV_reacspec/", celltype, "_seed_nodes_phuego.csv")
  output_file_txt <- paste0(directory_data,"/seed_nodes_scphuEGOTXT_reacspec/", celltype, "_seed_nodes_phuego.txt")
  if (!(file.exists(output_file_csv))){
    dir.create(file.path(directory_data,"seed_nodes_scphuEGOCSV"))
  }
  if (!(file.exists(output_file_txt))){
    dir.create(file.path(directory_data,"seed_nodes_scphuEGOTXT"))
  }
  
  saveRDS(merged_data, output_file_csv)
  write.table(merged_data[, c("Protein","RWRS")], output_file_txt, col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
  
  names(number_genetypes[[celltype]]) = c("DEGs","TFs","DErecs")
  return(number_genetypes)
}

# Loop over each cell type
type_genes <- List()
type_genes2 <- List()
DEGS <- c() #DEGs
TFS <- c() #TFs
DEREC <- c() #DE receptors
total <- c() #All genes selected
for (celltype in names(signature)) { #for each cell type
  types_thiscelltype <- process_cell_type(celltype)
  type_genes <- append(type_genes, types_thiscelltype)
  DEGS <- append(DEGS,type_genes[[celltype]][["DEGs"]])
  TFS <- append(TFS,type_genes[[celltype]][["TFs"]])
  DEREC <- append(DEREC,type_genes[[celltype]][["DErecs"]])
  total <- append(total, type_genes[[celltype]][["DEGs"]]+type_genes[[celltype]][["TFs"]]+type_genes[[celltype]][["DErecs"]])
}
type_genes2[["DEGS"]] <- DEGS
names(type_genes2[["DEGS"]]) <- names(signature)
type_genes2[["TFS"]] <- TFS
names(type_genes2[["TFS"]]) <- names(signature)
type_genes2[["DEREC"]] <- DEREC
names(type_genes2[["DEREC"]]) <- names(signature)
type_genes2[["total"]] <- total
names(type_genes2[["total"]]) <- names(signature)

hist(type_genes2$total,breaks=20,xlab="size of sc-phuEGO input seed nodes")

#Save the number of each gene type (DEG, DE receptor, TF) for each cell type
saveRDS(type_genes2,paste0(directory_output,"/genes_in_seednodes.csv"))
#Load the number of each gene type (DEG, DE receptor, TF) for each cell type
type_genes <- readRDS(paste0(directory_output,"/genes_in_seednodes.csv"))