
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5) {
  stop("Argument missing! Usage : scrip.R chrX start stop target_name path_to_file_list [path_to_snp_file] [highlight_region fac. start:end]")
}


current_chr = args[1]
current_start = as.numeric(args[2])
current_stop = as.numeric(args[3])
current_target = args[4]
file_list = args[5]
snp_file = args[6] # facultatif
highlight_region = args[7] # facultatif

suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbio))
#suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(biovizBase))
#suppressMessages(library(org.Hs.eg.db))
library(TxDb.Hsapiens.UCSC.hg19.knownGene, lib.loc="local_lib/" )
library(org.Hs.eg.db, lib.loc="local_lib/" )
library(parallel)

######################### FUNCTIONS ########################

parse_file_list = function(file_path) {
  
  fpe <- read.table(file_path, header=FALSE, stringsAsFactors=FALSE)
  l = list()
  n = length(fpe[[1]])
  for (i in(1:n)){
    exp = as.name(fpe[i,1])
    l <- c(l, list(strsplit(fpe[i,2], split = ",", fixed = TRUE)[[1]])) 
  }
  names(l)= fpe[[1]]
  
  l
}

file_path_sans_ext = function(bam_files) {
  files_name = c()
  for(i in (1 : length(bam_files))) {
    v=unlist(strsplit(bam_files[i],split="/",fixed=TRUE))
    m = sub(".bam","",v[length(v)])
    files_name = c(files_name,m)
  }
  files_name
}

######################### PROCESS ########################
border = 10000
current_range <- GRanges(current_chr, IRanges(current_start - border, current_stop + border))
bam_files_list <-  parse_file_list(file_list)

print(names(bam_files_list))

means =  list()

for (i in (1:length(bam_files_list))) {	
  bam_files_4_an_exp = bam_files_list[[i]]
  names(bam_files_4_an_exp) <- file_path_sans_ext(bam_files_4_an_exp)
  param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery=FALSE))
  counts <- mclapply(bam_files_4_an_exp, function(x) countBam(x, param = param)$records, mc.cores = 4)
  names(counts) <- names(bam_files_4_an_exp)
  param <- ScanBamParam(which = current_range)
  alignments <- mclapply(bam_files_4_an_exp, readGAlignments, param=param, mc.cores = 4)
  
  get_coverage <- function(bam_file) {
    weight <- 1 / (counts[[bam_file]] / 1000000)
    coverage <- coverage(alignments[[bam_file]], weight=weight)[current_range]
    coverage <- as.numeric(coverage[[1]])
    coverage[coverage < 0] <- 0
    coverage
  }
  
  coverages <- mclapply(names(bam_files_4_an_exp), get_coverage, mc.cores = 4)
  means <- c(means, list(colMeans(do.call("rbind", coverages))))
}

print("Means OK!")
if (is.na(highlight_region)) {
  highlight_region = "0:0"
}

names(means) = names(bam_files_list)
top_value <- max(unlist(mclapply(means, max, mc.cores = 4)))
trackname <- character()

get_plot <- function(bam_file_name) {
  coverage <- means[[bam_file_name]]
  position <- seq(start(current_range), end(current_range))
  trackname <- bam_file_name
  data <- data.frame(position = position, coverage = coverage)
  highlight_region = strsplit(x = highlight_region, fixed = T, split = ":")[[1]]
  d = data.frame(x1=as.numeric(highlight_region[1]), x2=as.numeric(highlight_region[2]), y1=0, y2=top_value)
  g = ggplot() + geom_line(data = data, mapping =  aes(x = position, y = coverage)) + geom_rect(data = d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha = 0.4) + ylim(0,top_value) + ylab(trackname) + theme_bw() + xlim(current_range)
  g
}

plots <- mclapply(names(means), get_plot, mc.cores = 4)

print("Plots OK!")

###################### SNPS DATA ###################### 
snps_track = NULL

if(!is.null(snp_file) && file.exists(snp_file)) {
  snps_df <- read.table(snp_file, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
  
  if(nrow(snps_df) > 0){
    ids <- as.character(snps_df$id)
    
    snps_ranges <- IRanges(as.numeric(as.character(snps_df$start)), as.numeric(as.character(snps_df$end)))
    snps <- GRanges(seqnames = as.character(seqnames(current_range)), ranges = snps_ranges, imp = snps_df$metadata)
    snps$name <- ids
    
    snps_track <- autoplot(snps, aes(color=imp)) +
      geom_text(aes(x = start, y = 1, label=name, angle = 90, vjust=-1), size = 1, color = "blue") +
      theme_bw() + xlim(current_range) + ylab("Variants") + guides(color= FALSE)
    
    print("Track snps -> DONE")
  }
}

###################### GENE DATA ###################### 
gr_txdb <- crunch(TxDb.Hsapiens.UCSC.hg19.knownGene, which = current_range)
colnames(values(gr_txdb))[4] <- "model"
gr_txdb$symbols <- select(org.Hs.eg.db,
                          keys = as.character(gr_txdb$gene_id),
                          column = "SYMBOL",
                          keytype = "ENTREZID")$SYMBOL
i <- which(gr_txdb$model == "gap"
)
gr_txdb <- gr_txdb[-i]
levels(gr_txdb) <- c("cds", "exon", "utr")
grl_txdb <- split(gr_txdb, gr_txdb$symbols)


genes <- autoplot(grl_txdb, aes(type = model, label = "Annotations")) + theme_bw() + xlim(current_range) + ylab("Annotations")

print("Track gene -> DONE")

if(is.null(snps_track)){
  plots = c(plots, genes)
} else {
  plots = c(plots, snps_track, genes)
}


pdf(paste0(current_target,"_rnaseq.pdf"))

tracks(plots, label.text.cex = 0.5) + (xlim(current_range - border))

dev.off()

save(means, file = paste0("means_",current_target,"_rnaseq.rda"))
save(plots, file = paste0("plots_",current_target,"_rnaseq_plots.rda"))


