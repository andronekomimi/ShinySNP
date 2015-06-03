############### LIBRARIES ###############

library(rhdf5)
library(parallel)

############### FUNCTIONS ###############
# transcript_id
df_col2transcript_id <- function(my.column) {
  my.str = unlist(strsplit(my.column, split = ";"))
  transcript_id_idx = grep(my.str, pattern = "transcript_id")
  transcript_id = gsub(pattern = " transcript_id ",x = my.str[transcript_id_idx], replacement = "")
  
  transcript_id
}

# transcript annotation
df_col2transcript_annot <- function(my.column) {
  my.str = unlist(strsplit(my.column, split = ";"))
  annotation_idx = grep(my.str, pattern = "tstatus")
  annotation = ! grepl(x = my.str[annotation_idx], pattern = "unannotated")
  
  annotation
}

############### PROCESS  ###############

# First : create the file
h5createFile("myhdf5file.h5")

# Second : create the file architecture
# local/cell/chr/experiments/OBJECT
# public/cell/chr/experiments/OBJECT
# common/chr/experiments/OBJECT
# common/misc/OBJECT

my.datasets = c("interne","externe")
cells = c("HMEC","K562","MCF7")
# chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
#          "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
#          "chr18","chr19","chr20","chr21","chr22", "chrX", "chrY", "chrM")
chrs = c("chr6","chr7")

names(chrs) = chrs 

methods = c("3C","4C","5C","ChIA-PET","Hi-C","IM-PET")
common_methods = c("ENHANCER", "SUPER-ENHANCER", "LNCRNA")

names(methods) = methods

for (dataset in my.datasets) {
  
  h5createGroup("myhdf5file.h5",dataset)
  
  for (cell in cells) {
    
    h5createGroup("myhdf5file.h5",paste0(dataset,"/",cell))
    
    for (chr in chrs) {
      
      h5createGroup("myhdf5file.h5",paste0(dataset,"/",cell,"/",chr))
      
      for (method in methods) {
        h5createGroup("myhdf5file.h5",paste0(dataset,"/",cell,"/",chr,"/",method))
      }      
    }
  }
}

h5createGroup("myhdf5file.h5","common")
h5createGroup("myhdf5file.h5",paste0("common/misc"))

for (chr in chrs) {
  h5createGroup("myhdf5file.h5",paste0("common/",chr))
  
  for (method in common_methods) {
    
    h5createGroup("myhdf5file.h5",paste0("common/",chr,"/",method))
  }
  
}



h5ls("myhdf5file.h5")

# Third : remplissage avec data subset interne specific
############################ CHIAPET ############################ 
############################ SINGAPORE DATA ############################  
chiapet_k562_lane13_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane13_"
chiapet_k562_lane24_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane24_"
chiapet_mcf7_lane11_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane11_"
chiapet_mcf7_lane23_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane23_"
############################ STANDFORD DATA ############################  
chiapet_k562_rep1_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep1_"
chiapet_k562_rep2_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep2_"
chiapet_mcf7_rep3_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep3_"
chiapet_mcf7_rep4_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep4_"
############################ HIC ############################ 
hic_hmec_file_1 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_HMEC_Arrowhead_domainlist_"
hic_hmec_file_2 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_HMEC_HiCCUPS_looplist_"
hic_k562_file_1 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_K562_Arrowhead_domainlist_"
hic_k562_file_2 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_K562_HiCCUPS_looplist_"


sgp_k562_files_path = c(chiapet_k562_lane13_path, chiapet_k562_lane24_path)

sgp_mcf7_files_path = c(chiapet_mcf7_lane11_path, chiapet_mcf7_lane23_path)

sdf_k562_file_path = c(chiapet_k562_rep1_path, chiapet_k562_rep2_path)

sdf_mcf7_file_path = c(chiapet_mcf7_rep3_path, chiapet_mcf7_rep4_path)

hic_hmec_files_path = c(hic_hmec_file_1, hic_hmec_file_2)

hic_k562_files_path = c(hic_k562_file_1, hic_k562_file_2)

print("start interne data loading")

for (chr in chrs) {
  #Chiapet
  #K562
  #file1 de sgp
  
  my.file <- paste0(sgp_k562_files_path[1], chr)
  df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df1 <- subset(df1, select=c(1,2,3,4,5,6,7))
  colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  #file2 de sgp
  my.file <- paste0(sgp_k562_files_path[2], chr)
  df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df2 <- subset(df2, select=c(1,2,3,4,5,6,7))
  colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  #file1 de sfd
  my.file <- paste0(sdf_k562_file_path[1], chr)
  df3 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  my.data <- strsplit(gsub("\\.\\.", "-", df3[,4]), '[:,-]')
  df3 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE))
  colnames(df3) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  #file2 de sfd
  my.file <- paste0(sdf_k562_file_path[2], chr)
  df4 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  my.data <- strsplit(gsub("\\.\\.", "-", df4[,4]), '[:,-]')
  df4 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE))
  colnames(df4) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  # DF containing all chiapet data for k562
  k562 <- rbind(df1, df2, df3, df4)
  
  if(nrow(k562) > 0) {
    k562 <- transform(k562, InteractorAStart = as.numeric(InteractorAStart), 
                      InteractorAEnd = as.numeric(InteractorAStart),
                      InteractorBStart = as.numeric(InteractorAEnd),
                      InteractorBEnd = as.numeric(InteractorBEnd))
    k562 <- k562[ order(k562$InteractorAStart),]
    
    path <- paste0("interne/K562/",chr,"/ChIA-PET/k562")
    h5write(k562, "myhdf5file.h5",path)
  }
  
  #MCF7
  #file1 de sgp
  my.file <- paste0(sgp_mcf7_files_path[1], chr)
  df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df1 <- subset(df1, select=c(1,2,3,4,5,6,7))
  colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  #file2 de sgp
  my.file <- paste0(sgp_mcf7_files_path[2], chr)
  df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df2 <- subset(df2, select=c(1,2,3,4,5,6,7))
  colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  #file1 de sfd
  my.file <- paste0(sdf_mcf7_file_path[1], chr)
  df3 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  my.data <- strsplit(gsub("\\.\\.", "-", df3[,4]), '[:,-]')
  df3 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE))
  colnames(df3) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  #file2 de sfd
  my.file <- paste0(sdf_mcf7_file_path[2], chr)
  df4 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  my.data <- strsplit(gsub("\\.\\.", "-", df4[,4]), '[:,-]')
  df4 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE))
  colnames(df4) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  # DF containing all chiapet data for mcf7
  mcf7 <- rbind(df1, df2, df3, df4)
  if(nrow(mcf7) > 0) {
    mcf7 <- transform(mcf7, InteractorAStart = as.numeric(InteractorAStart), 
                      InteractorAEnd = as.numeric(InteractorAStart),
                      InteractorBStart = as.numeric(InteractorAEnd),
                      InteractorBEnd = as.numeric(InteractorBEnd))
    mcf7 <- mcf7[ order(mcf7$InteractorAStart),]
    
    path <- paste0("interne/MCF7/",chr,"/ChIA-PET/mcf7")
    h5write(mcf7, "myhdf5file.h5",path)
  }
  
  #HIC
  #HMEC
  #file1 de Arrow
  my.file <- paste0(hic_hmec_files_path[1], chr)
  df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df1 <- subset(df1, select=c(1,2,3,4,5,6,8))
  colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  #file2 de cups
  my.file <- paste0(hic_hmec_files_path[2], chr)
  df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df2 <- subset(df2, select=c(1,2,3,4,5,6,8))
  colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  # DF containing all hic data for hmec
  hmec <- rbind(df1, df2)
  if(nrow(hmec) > 0) {
    
    hmec <- transform(hmec, InteractorAStart = as.numeric(InteractorAStart), 
                      InteractorAEnd = as.numeric(InteractorAStart),
                      InteractorBStart = as.numeric(InteractorAEnd),
                      InteractorBEnd = as.numeric(InteractorBEnd))
    
    hmec <- hmec[ order(hmec$InteractorAStart),]
    
    path <- paste0("interne/HMEC/",chr,"/Hi-C/hmec")
    h5write(hmec, "myhdf5file.h5",path)
  }
  
  #K562
  #file1 de Arrow
  my.file <- paste0(hic_k562_files_path[1], chr)
  df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df1 <- subset(df1, select=c(1,2,3,4,5,6,8))
  colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  #file2 de cups
  my.file <- paste0(hic_k562_files_path[2], chr)
  df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
  df2 <- subset(df2, select=c(1,2,3,4,5,6,8))
  colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                     "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  
  # DF containing all hic data for hmec
  k562 <- rbind(df1, df2)
  if(nrow(k562) > 0) {
    k562 <- transform(k562, InteractorAStart = as.numeric(InteractorAStart), 
                      InteractorAEnd = as.numeric(InteractorAStart),
                      InteractorBStart = as.numeric(InteractorAEnd),
                      InteractorBEnd = as.numeric(InteractorBEnd))
    k562 <- k562[ order(k562$InteractorAStart),]
    
    path <- paste0("interne/K562/",chr,"/Hi-C/k562")
    h5write(hmec, "myhdf5file.h5",path)
  }
  
}

# Fourth : remplissage avec data subset externe specific
# 4DGENOME
print("start externe data loading")
data_file <- "/home/nekomimi/Workspace/SNPVIZU/data/4DGenome_HomoSapiens_hg19.txt"
df1 <- read.table(data_file, header=TRUE, stringsAsFactors=TRUE, sep = "\t")

for (chr in chrs) {
  chr_df <- subset(df1, df1$InteractorAChr == chrs[[chr]])
  chr_df = chr_df[with(chr_df, order(chr_df$InteractorAStart)),]
  
  for(method in methods) {
    temp = subset(chr_df, chr_df$Detection_Method == method)
    
    mcf7 = subset(temp, temp$Cell.Tissue == "MCF7")
    k562 = subset(temp, temp$Cell.Tissue == "K562")
    hmec = subset(temp, temp$Cell.Tissue == "HMEC")
    
    if(nrow(mcf7) > 0) {
      mcf7 <- transform(mcf7, InteractorAStart = as.numeric(InteractorAStart), 
                        InteractorAEnd = as.numeric(InteractorAStart),
                        InteractorBStart = as.numeric(InteractorAEnd),
                        InteractorBEnd = as.numeric(InteractorBEnd))
      mcf7 <- mcf7[ order(mcf7$InteractorAStart),]
      
      path <- paste0("externe/MCF7/",chr,"/",method,"/mcf7")
      h5write(mcf7, "myhdf5file.h5",path)
    }
    
    if(nrow(k562) > 0) {
      k562 <- transform(k562, InteractorAStart = as.numeric(InteractorAStart), 
                        InteractorAEnd = as.numeric(InteractorAStart),
                        InteractorBStart = as.numeric(InteractorAEnd),
                        InteractorBEnd = as.numeric(InteractorBEnd))
      k562 <- k562[ order(k562$InteractorAStart),]
      
      path <- paste0("externe/K562/",chr,"/",method,"/k562")
      h5write(k562, "myhdf5file.h5",path)
    }
    
    if(nrow(hmec) > 0) {
      hmec <- transform(hmec, InteractorAStart = as.numeric(InteractorAStart), 
                        InteractorAEnd = as.numeric(InteractorAStart),
                        InteractorBStart = as.numeric(InteractorAEnd),
                        InteractorBEnd = as.numeric(InteractorBEnd))
      hmec <- hmec[ order(hmec$InteractorAStart),]
      
      path <- paste0("externe/HMEC/",chr,"/",method,"/hmec")
      h5write(hmec, "myhdf5file.h5",path)
    }
  }
}

# Fifth : remplissage avec data communes
############################ LNCRNA ############################ 

print("start common data loading")

files_path <- "/home/nekomimi/Workspace/COLLAB/mitranscriptome.gtf/mitranscriptome_"

for (chr in chrs) {
  my.file <- paste0(files_path, chr)
  lncrna <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
  colnames(lncrna) = c("InteractorAChr","source","type", "InteractorAStart",
                   "InteractorAEnd","metric","strand", "transcript_id",
                   "annotation")
  lncrna <- subset(lncrna,lncrna$type == "transcript")
  lncrna <- transform(lncrna, transcript_id = df_col2transcript_id(lncrna$annotation))
  lncrna <- transform(lncrna, annotation = df_col2transcript_id(lncrna$annotation))
  lncrna <- transform(lncrna, InteractorAStart = as.numeric(InteractorAStart), 
                  InteractorAEnd = as.numeric(InteractorAStart), 
                  annotation = as.logical(annotation))
  path <- paste0("common/",chr,"/LNCRNA/lncrna")
  h5write(lncrna, "myhdf5file.h5",path)
}


