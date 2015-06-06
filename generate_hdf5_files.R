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
t0 = Sys.time()

# First : create the file
h5createFile("data/data.h5")

# Second : create the file architecture
# local/cell/chr/experiments/OBJECT
# public/cell/chr/experiments/OBJECT
# common/cell/chr/experiement/OBJECT
# common/chr/experiments/OBJECT
# common/misc/experiment/OBJECT

my.datasets = c("interne","externe", "common")

cells = c("HMEC","K562","MCF7")
chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
         "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
         "chr18","chr19","chr20","chr21","chr22", "chrX", "chrY", "chrM")
# chrs = c("chr6","chr7")

names(chrs) = chrs 

methods = c("3C","4C","5C","ChIA-PET","Hi-C","IM-PET")
common_methods = c("ENHANCER", "LNCRNA")
chiapet_subcat = c("ER","CTCF")

names(methods) = methods

for (dataset in my.datasets) {
  
  h5createGroup("data/data.h5",dataset)
  
  for (cell in cells) {
    
    h5createGroup("data/data.h5",paste0(dataset,"/",cell))
    
    for (chr in chrs) {
      
      h5createGroup("data/data.h5",paste0(dataset,"/",cell,"/",chr))
      
      for (method in methods) {
        h5createGroup("data/data.h5",paste0(dataset,"/",cell,"/",chr,"/",method))
      }
      
      for (sub in chiapet_subcat) {
        h5createGroup("data/data.h5",paste0(dataset,"/",cell,"/",chr,"/ChIA-PET/",sub))
      }
      
      h5createGroup("data/data.h5",paste0(dataset,"/",cell,"/",chr,"/SUPER-ENHANCER"))
      
    }
  }
}

h5createGroup("data/data.h5",paste0("common/misc"))

for (chr in chrs) {
  h5createGroup("data/data.h5",paste0("common/",chr))
  
  for (method in common_methods) {
    
    h5createGroup("data/data.h5",paste0("common/",chr,"/",method))
  }
  
}



h5ls("data/data.h5")

# Third : remplissage avec data subset interne specific
############################ CHIAPET ############################ 
############################ SINGAPORE DATA ############################  
chiapet_k562_lane13_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane13_"
chiapet_k562_lane24_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane24_"
chiapet_mcf7_lane11_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane11_"
chiapet_mcf7_lane23_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane23_"
chiapet_mcf7_er_file1 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_ER_rep1_"
chiapet_mcf7_er_file2 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_ER_rep2_"
#chiapet_mcf7_er_file3 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_ER_rep3_"
chiapet_mcf7_ctcf_file1 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_CTCF_rep1_"
chiapet_mcf7_ctcf_file2 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_CTCF_rep2_"
chiapet_k562_ctcf_file1 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_CTCF_rep1_"
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

sgp_mcf7_er_files_path = c(chiapet_mcf7_er_file1, chiapet_mcf7_er_file2)

sgp_mcf7_ctcf_files_path = c(chiapet_mcf7_ctcf_file1, chiapet_mcf7_ctcf_file2)

sgp_k562_ctcf_files_path = c(chiapet_k562_ctcf_file1)

sdf_k562_file_path = c(chiapet_k562_rep1_path, chiapet_k562_rep2_path)

sdf_mcf7_file_path = c(chiapet_mcf7_rep3_path, chiapet_mcf7_rep4_path)

hic_hmec_files_path = c(hic_hmec_file_1, hic_hmec_file_2)

hic_k562_files_path = c(hic_k562_file_1, hic_k562_file_2)

print("start interne data loading")

empty_df1 <- data.frame(InteractorAChr= character(),
                        InteractorAStart=integer(),
                        InteractorAEnd=integer(),
                        InteractorBChr= character(),
                        InteractorBStart=integer(),
                        InteractorBEnd=integer(),
                        metric=character(), 
                        stringsAsFactors=FALSE)

empty_df2 <- data.frame(InteractorAChr= character(),
                        InteractorAStart=integer(),
                        InteractorBEnd=integer(),
                        peak= character(),
                        metric=character(), 
                        stringsAsFactors=FALSE)


for (chr in chrs) {
  #Chiapet
  #K562
  #file1 de sgp
  my.file <- paste0(sgp_k562_files_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df1 <- subset(df1, select=c(1,2,3,4,5,6,7))
    colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('3.EMPTY')
    df1 <- empty_df1
  }
  
  #file2 de sgp
  my.file <- paste0(sgp_k562_files_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df2 <- subset(df2, select=c(1,2,3,4,5,6,7))
    colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('4.EMPTY')
    df2 <- empty_df1
  }
  
  #file1 de sfd
  my.file <- paste0(sdf_k562_file_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df3 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df3[,4]), '[:,-]')
    df3 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df3) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('5.EMPTY')
    df3 <- empty_df1
  }
  
  #file2 de sfd
  my.file <- paste0(sdf_k562_file_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df4 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df4[,4]), '[:,-]')
    df4 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df4) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('6.EMPTY')
    df4 <- empty_df1
  }
  
  # DF containing all chiapet data for k562
  #   print(names(df1))
  #   print(head(df2))
  #   print(names(df3))
  #   print(names(df4))
  #   
  df <- rbind(df1, df2, df3, df4)
  
  if(nrow(df) > 0) {
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart),
                    InteractorBStart = as.numeric(InteractorAEnd),
                    InteractorBEnd = as.numeric(InteractorBEnd))
    df <- df[ order(df$InteractorAStart),]
    
    path <- paste0("interne/K562/",chr,"/ChIA-PET/df")
    h5write(df, "data/data.h5",path)
  }
  
  #MCF7
  #file1 de sgp
  my.file <- paste0(sgp_mcf7_files_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df1 <- subset(df1, select=c(1,2,3,4,5,6,7))
    colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('7.EMPTY')
    df1 <- empty_df1
  }
  
  
  #file2 de sgp
  my.file <- paste0(sgp_mcf7_files_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df2 <- subset(df2, select=c(1,2,3,4,5,6,7))
    colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('8.EMPTY')
    df2 <- empty_df1
  }
  
  #file1 de sfd
  my.file <- paste0(sdf_mcf7_file_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df3 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df3[,4]), '[:,-]')
    df3 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df3) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('9.EMPTY')
    df3 <- empty_df1
  }
  
  
  
  #file2 de sfd
  my.file <- paste0(sdf_mcf7_file_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df4 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df4[,4]), '[:,-]')
    df4 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df4) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
    
  }
  else
  {
    #print('10.EMPTY')
    df4 <- empty_df1
  }
  
  
  # DF containing all chiapet data for mcf7
  df <- rbind(df1, df2, df3, df4)
  if(nrow(df) > 0) {
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart),
                    InteractorBStart = as.numeric(InteractorAEnd),
                    InteractorBEnd = as.numeric(InteractorBEnd))
    df <- df[ order(df$InteractorAStart),]
    
    path <- paste0("interne/MCF7/",chr,"/ChIA-PET/df")
    h5write(df, "data/data.h5",path)
  }
  
  #HIC
  #HMEC
  #file1 de Arrow
  my.file <- paste0(hic_hmec_files_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df1 <- subset(df1, select=c(1,2,3,4,5,6,8))
    colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('11.EMPTY')
    df1 <- empty_df1
  }
  
  
  #file2 de cups
  my.file <- paste0(hic_hmec_files_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df2 <- subset(df2, select=c(1,2,3,4,5,6,8))
    colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('12.EMPTY')
    df2 <- empty_df1
  }
  
  
  # DF containing all hic data for hmec
  df <- rbind(df1, df2)
  if(nrow(df) > 0) {
    
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart),
                    InteractorBStart = as.numeric(InteractorAEnd),
                    InteractorBEnd = as.numeric(InteractorBEnd))
    
    df <- df[ order(df$InteractorAStart),]
    
    path <- paste0("interne/HMEC/",chr,"/Hi-C/df")
    h5write(df, "data/data.h5",path)
  }
  
  #K562
  #file1 de Arrow
  my.file <- paste0(hic_k562_files_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df1 <- subset(df1, select=c(1,2,3,4,5,6,8))
    colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('1.EMPTY')
    df1 <- empty_df1
  }
  
  
  
  #file2 de cups
  my.file <- paste0(hic_k562_files_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    df2 <- subset(df2, select=c(1,2,3,4,5,6,8))
    colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('2.EMPTY')
    df2 <- empty_df1
  }
  
  
  # DF containing all hic data for hmec
  df <- rbind(df1, df2)
  if(nrow(df) > 0) {
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart),
                    InteractorBStart = as.numeric(InteractorAEnd),
                    InteractorBEnd = as.numeric(InteractorBEnd))
    df <- df[ order(df$InteractorAStart),]
    
    path <- paste0("interne/K562/",chr,"/Hi-C/df")
    h5write(df, "data/data.h5",path)
  }
  
  # CHIAPET ER MCF7
  my.file <- paste0(sgp_mcf7_er_files_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df1[,4]), '[:,-]')
    df1 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('3.EMPTY')
    df1 <- empty_df1
  }
  
  my.file <- paste0(sgp_mcf7_er_files_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df2[,4]), '[:,-]')
    df2 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('3.EMPTY')
    df2 <- empty_df1
  }
  
#   my.file <- paste0(sgp_mcf7_er_files_path[3], chr)
#   #print(my.file)
#   if(file.info(my.file)$size > 0) {
#     df3 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
#     my.data <- strsplit(gsub("\\.\\.", "-", df2[,4]), '[:,-]')
#     df3 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
#     colnames(df3) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
#                        "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
#   }
#   else
#   {
#     #print('3.EMPTY')
#     df3 <- empty_df1
#   }  
  
  df <- rbind(df1, df2)
  
  #print(head(df))
  
  if(nrow(df) > 0) {
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart),
                    InteractorBStart = as.numeric(InteractorAEnd),
                    InteractorBEnd = as.numeric(InteractorBEnd))
    df <- df[ order(df$InteractorAStart),]
    
    path <- paste0("interne/MCF7/",chr,"/ChIA-PET/ER/df")
    h5write(df, "data/data.h5",path)
  }
  
  # CHIAPET CTCF MCF7
  my.file <- paste0(sgp_mcf7_ctcf_files_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df1 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df1[,4]), '[:,-]')
    df1 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df1) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('3.EMPTY')
    df1 <- empty_df1
  }
  
  my.file <- paste0(sgp_mcf7_ctcf_files_path[2], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df2 <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df2[,4]), '[:,-]')
    df2 =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df2) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('3.EMPTY')
    df2 <- empty_df1
  }
  
  df <- rbind(df1, df2)
  
  #print(head(df))
  
  if(nrow(df) > 0) {
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart),
                    InteractorBStart = as.numeric(InteractorAEnd),
                    InteractorBEnd = as.numeric(InteractorBEnd))
    df <- df[ order(df$InteractorAStart),]
    
    path <- paste0("interne/MCF7/",chr,"/ChIA-PET/CTCF/df")
    h5write(df, "data/data.h5",path)
  }
  
  # CHIAPET CTCF K562
  my.file <- paste0(sgp_k562_ctcf_files_path[1], chr)
  #print(my.file)
  if(file.info(my.file)$size > 0) {
    df <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    my.data <- strsplit(gsub("\\.\\.", "-", df[,4]), '[:,-]')
    df =  data.frame(matrix(unlist(my.data), ncol = 7, byrow = TRUE), stringsAsFactors=FALSE)
    colnames(df) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                       "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
  }
  else
  {
    #print('3.EMPTY')
    df <- empty_df1
  }
  
  #print(head(df))
  
  if(nrow(df) > 0) {
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart),
                    InteractorBStart = as.numeric(InteractorAEnd),
                    InteractorBEnd = as.numeric(InteractorBEnd))
    df <- df[ order(df$InteractorAStart),]
    
    path <- paste0("interne/K562/",chr,"/ChIA-PET/CTCF/df")
    h5write(df, "data/data.h5",path)
  }
}

# Fourth : remplissage avec data subset externe specific
# 4DGENOME
print("start externe data loading")
data_file <- "/home/nekomimi/Workspace/SNPVIZU/data/4DGenome_HomoSapiens_hg19.txt"
df1 <- read.table(data_file, header=TRUE, stringsAsFactors=FALSE, sep = "\t")

for (chr in chrs) {
  chr_df <- subset(df1, df1$InteractorAChr == chrs[[chr]])
  chr_df = chr_df[with(chr_df, order(chr_df$InteractorAStart)),]
  
  for(method in methods) {
    temp = subset(chr_df, chr_df$Detection_Method == method)
    
    df = subset(temp, temp$Cell.Tissue == "MCF7")
    if(nrow(df) > 0) {
      df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                      InteractorAEnd = as.numeric(InteractorAStart),
                      InteractorBStart = as.numeric(InteractorAEnd),
                      InteractorBEnd = as.numeric(InteractorBEnd))
      df <- df[ order(df$InteractorAStart),]
      
      path <- paste0("externe/MCF7/",chr,"/",method,"/df")
      h5write(df, "data/data.h5",path)
    }
    
    df = subset(temp, temp$Cell.Tissue == "K562")
    if(nrow(df) > 0) {
      df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                      InteractorAEnd = as.numeric(InteractorAStart),
                      InteractorBStart = as.numeric(InteractorAEnd),
                      InteractorBEnd = as.numeric(InteractorBEnd))
      df <- df[ order(df$InteractorAStart),]
      
      path <- paste0("externe/K562/",chr,"/",method,"/df")
      h5write(df, "data/data.h5",path)
    }
    
    df = subset(temp, temp$Cell.Tissue == "HMEC")
    if(nrow(df) > 0) {
      df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                      InteractorAEnd = as.numeric(InteractorAStart),
                      InteractorBStart = as.numeric(InteractorAEnd),
                      InteractorBEnd = as.numeric(InteractorBEnd))
      df <- df[ order(df$InteractorAStart),]
      
      path <- paste0("externe/HMEC/",chr,"/",method,"/df")
      h5write(df, "data/data.h5",path)
    }
  }
}

# Fifth : remplissage avec data communes
############################ LNCRNA ############################ 

print("start common data loading")

files_path <- "/home/nekomimi/Workspace/COLLAB/mitranscriptome.gtf/mitranscriptome_"

for (chr in chrs) {
  my.file <- paste0(files_path, chr)
  df <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
  if(nrow(df) > 0){
    colnames(df) = c("InteractorAChr","source","type", "InteractorAStart",
                     "InteractorBEnd","metric","strand", "transcript_id",
                     "annotation")
    df <- subset(df,df$type == "transcript")
    df <- transform(df, transcript_id = df_col2transcript_id(df$annotation))
    df <- transform(df, annotation = df_col2transcript_id(df$annotation))
    df <- transform(df, InteractorAStart = as.numeric(InteractorAStart), 
                    InteractorAEnd = as.numeric(InteractorAStart), 
                    annotation = as.logical(annotation))
    path <- paste0("common/",chr,"/LNCRNA/df")
    h5write(df, "data/data.h5",path)
  }
}

lncrna_expr_file <- "/home/nekomimi/Workspace/COLLAB/mitranscriptome.expr.fpkm_select.tsv"
df <- read.table(lncrna_expr_file, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
h5createGroup("data/data.h5",paste0("common/misc/LNCRNA-EXP"))
h5write(df, "data/data.h5","common/misc/LNCRNA-EXP/df")

############################ ENHANCERS & SUPER ENHANCERS ############################ 

enhancers_file <- "/home/nekomimi/Workspace/SNPVIZU/data/ENHANCERS/human_permissive_enhancers_phase_1_and_2.bed_"
super_enhancers_file <- "/home/nekomimi/Workspace/SNPVIZU/data/ENHANCERS/super_enhancers_"

for(chr in chrs) {
  my.file <- paste0(enhancers_file, chr)
  #print(my.file)
  if(file.info(my.file)$size > 0){
    df = read.csv(my.file, sep="\t", header= FALSE, 
                  stringsAsFactors=FALSE)
    df = subset(df, select=c(1,2,3,1,7,8,5))
    colnames(df) <- c("InteractorAChr","InteractorAStart","InteractorAEnd",
                      "InteractorBChr","InteractorBStart","InteractorBEnd","metric")
    if(nrow(df) > 0 ){ 
      path <- paste0("common/",chr,"/ENHANCER/df")
      h5write(df, "data/data.h5",path)
    }
  }
  
  my.file = paste0(super_enhancers_file,"MCF7.bed_", chr)
  #print(my.file)
  if(file.info(my.file)$size > 0){
    df = read.csv(my.file, sep="\t", header= FALSE, 
                  stringsAsFactors=FALSE)
    
    if(nrow(df) > 0 ){ 
      colnames(df) <- c("InteractorAChr","InteractorAStart","InteractorBEnd",
                        "peak","metric")
      path <- paste0("common/MCF7/",chr,"/SUPER-ENHANCER/df")
      h5write(df, "data/data.h5",path)
    }
  }
  
  my.file = paste0(super_enhancers_file,"K562.bed_", chr)
  #print(my.file)
  if(file.info(my.file)$size > 0){
    df = read.csv(my.file, sep="\t", header= FALSE, 
                  stringsAsFactors=FALSE)
    
    if(nrow(df) > 0 ){
      colnames(df) <- c("InteractorAChr","InteractorAStart","InteractorBEnd",
                        "peak","metric")
      path <- paste0("common/K562/",chr,"/SUPER-ENHANCER/df")
      h5write(df, "data/data.h5",path)
    }
  }
  
  my.file = paste0(super_enhancers_file,"HMEC.bed_", chr)
  #print(my.file)
  if(file.info(my.file)$size > 0){
    df = read.csv(my.file, sep="\t", header= FALSE, 
                  stringsAsFactors=FALSE)
    if(nrow(df) > 0){
      colnames(df) <- c("InteractorAChr","InteractorAStart","InteractorBEnd",
                        "peak","metric")
      path <- paste0("common/HMEC/",chr,"/SUPER-ENHANCER/df")
      h5write(df, "data/data.h5",path)
    }
  }
}

print(Sys.time()-t0)



