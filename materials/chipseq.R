t0 = Sys.time()
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5) {
  stop("Argument missing! Usage : scrip.R chrX start stop uniq_id path_to_file_list [figure_title] [path_to_snp_file] [path_to_highlight_file]")
}


current_chr = args[1]
current_start = as.numeric(args[2])
current_stop = as.numeric(args[3])
uniq_id = args[4]
file_list = args[5]
fig_title = args[6] # facultatif
snp_file = args[7] # facultatif
highlight_file = args[8] # facultatif

if(is.na(fig_title)) {
  fig_title = "RNA-Seq figure"
}

suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbio))
#suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
library(TxDb.Hsapiens.UCSC.hg19.knownGene, lib.loc="local_lib/" )
suppressMessages(library(biovizBase))
#suppressMessages(library(org.Hs.eg.db))
library(org.Hs.eg.db, lib.loc="local_lib/" )
library(parallel)

######################### FUNCTIONS ########################

parse_file_list = function(file_path) {
  
  fpe <- read.table(file_path, header=FALSE, stringsAsFactors=FALSE)
  l = list()
  ctr = c()
  n = length(fpe[[1]])
  for (i in(1:n)){
    l <- c(l, list(strsplit(fpe[i,2], split = ",", fixed = TRUE)[[1]]))
    ctr <- c(ctr,fpe[i,3])
  }
  names(l)= fpe[[1]]
  names(ctr)= fpe[[1]]
  
  list(l,ctr)
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

create_plots_list = function(myplots) {
  
  mylist = list()
  
  #EZH2
  if(! is.null(myplots$EZH2)) {
    myplot = myplots$EZH2
    myplot$layers[[1]]$geom_params$colour = "gray"
    mylist = c(mylist, list(myplot))
  }
  
  #CTCF
  if(! is.null(myplots$CTCF)) {
    myplot = myplots$CTCF
    myplot$layers[[1]]$geom_params$colour = "blue"
    mylist = c(mylist, list(myplot))
  }
  
  #H3K27me3
  if(! is.null(myplots$H3K27me3)) {
    myplot = myplots$H3K27me3
    myplot$layers[[1]]$geom_params$colour = "purple"
    myplot$scales$scales[[1]]$limits[2] = 1.2
    mylist = c(mylist, list(myplot))
  }
  
  #H3K4me1
  if(! is.null(myplots$H3K4me1)) {
    myplot = myplots$H3K4me1
    myplot$layers[[1]]$geom_params$colour = "#FAD400"
    myplot$scales$scales[[1]]$limits[2] = 1.2
    mylist = c(mylist, list(myplot))
  }
  
  
  #H3K27ac
  if(! is.null(myplots$H3K27ac)) {
    myplot = myplots$H3K27ac
    myplot$layers[[1]]$geom_params$colour = "#FAC000"
    myplot$scales$scales[[1]]$limits[2] = 1.2
    mylist = c(mylist, list(myplot))
  }
  
  #H3K4me2
  if(! is.null(myplots$H3K4me2)) {
    myplot = myplots$H3K4me2
    myplot$layers[[1]]$geom_params$colour = "#FA7D00"
    myplot$scales$scales[[1]]$limits[2] = 1.2
    mylist = c(mylist, list(myplot))
  }
  
  #H3K4me3
  if(! is.null(myplots$H3K4me3)) {
    myplot = myplots$H3K4me3
    myplot$layers[[1]]$geom_params$colour = "#FA7D00"
    myplot$scales$scales[[1]]$limits[2] = 1.2
    mylist = c(mylist, list(myplot))
  }
  
  
  #H3K36me3
  if(! is.null(myplots$H3K36me3)) {
    myplot = myplots$H3K36me3
    myplot$layers[[1]]$geom_params$colour = "#19C910"
    myplot$scales$scales[[1]]$limits[2] = 1.2
    mylist = c(mylist, list(myplot))
  }
  
  #rnaseq
  if(! is.null(myplots$rnaseq)) {
    myplot = myplots$rnaseq
    myplot$layers[[1]]$geom_params$colour = "#37A331"
    mylist = c(mylist, list(myplot))
  }
  
  if(! is.null(myplots$snps)) {
    myplot = myplots$snps
    mylist = c(mylist, list(myplot))
  }
  
  if(! is.null(myplots$annotation)) {
    myplot = myplots$annotation
    mylist = c(mylist, list(myplot))
  }
  
  mylist
}


######################### PROCESS ########################

current_range <- GRanges(current_chr, IRanges(current_start, current_stop))

ret <- parse_file_list(file_list)
bam_files_list <- ret[[1]]
controls <- ret[[2]]


print(paste0("bam_files_list length : ",length(bam_files_list)))

means = list()

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

top_value <- max(unlist(mclapply(means, max, mc.cores = 4)))
names(means) = names(bam_files_list)

roll_mean <- function(experiment, n=100){
  res = means[[experiment]]
  for(i in n:length(means[[experiment]])){
    res[i] = mean( means[[experiment]][(i-n):i])
  }
  res
}

apply_control <- function(experiment) {
  ctr =  controls[[experiment]]
  if (! is.na(ctr)) {
    new_mean = means[[experiment]] - means[[ctr]]
  }
  else
  {
    new_mean = means[[experiment]]
  }
  
  new_mean
}

new_means = mclapply(names(means), apply_control, mc.cores = 4)
new_means = mclapply(names(means), roll_mean, mc.cores = 4)
names(new_means) = names(means)
means = new_means

print(names(means))

trackname <- character()

get_plot <- function(bam_file_name) {
  coverage <- means[[bam_file_name]]
  position <- seq(start(current_range), end(current_range))
  trackname <- bam_file_name
  data <- data.frame(position = position, coverage = coverage)
  
  if(!is.null(highlight_file) && file.exists(highlight_file)) {
    hgs_df <- read.table(highlight_file, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
    #d <- data.frame(x1=as.numeric(hgs_df$start), x2=as.numeric(hgs_df$end), y1=0, y2=max(data$coverage), col = hgs_df$color)
    d <- data.frame(x1=as.numeric(hgs_df$start), x2=as.numeric(hgs_df$end), y1=0, y2=1, col = hgs_df$color)
    my.colors = d$col
    names(my.colors) = my.colors
  } else {
    d <- data.frame(x1=c(0), x2=c(0), y1=0, y2=top_value, col=c("black"))
  }
  
  g = ggplot() + geom_line(data = data, mapping =  aes(x = position, y = coverage)) + 
    geom_rect(data = d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.4, fill = my.colors) +
    ylim(0,top_value) + ylab(trackname) + theme_bw() + 
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90), axis.text.y = element_text(size = rel(0.5)))+ 
    xlim(current_range)
  
  print(paste0("Track plot for ", trackname," -> DONE"))
  
  g
}

plots <- mclapply(names(means), get_plot, mc.cores = 4)
names(plots) = names(means)


###################### SNPS DATA ###################### 
snps_track = NULL

if(!is.null(snp_file) && file.exists(snp_file)) {
  snps_df <- read.table(snp_file, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
  
  if(nrow(snps_df) > 0){
    ids <- as.character(snps_df$id)
    
    snps_ranges <- IRanges(as.numeric(as.character(snps_df$start)), as.numeric(as.character(snps_df$end)))
    snps <- GRanges(seqnames = as.character(seqnames(current_range)), ranges = snps_ranges, imp = snps_df$metadata)
    snps$name <- ids
    
    snps_track <- autoplot(snps, aes(fill=imp)) +
      geom_text(aes(x = start, y = 1, label=name, angle = 90, vjust=-1), size = 1, color = "blue") +
      theme_bw() + xlim(current_range) + ylab("Variants") + guides(fill= FALSE)
    
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
i <- which(gr_txdb$model == "gap")
gr_txdb <- gr_txdb[-i]
levels(gr_txdb) <- c("cds", "exon", "utr")
grl_txdb <- split(gr_txdb, gr_txdb$symbols)


genes <- autoplot(grl_txdb, aes(type = model, label = "Annotations")) + theme_bw() + xlim(current_range) + ylab("Annotations") + theme(axis.title.y = element_text(size = rel(0.5), angle = 90), axis.text.y = element_text(size = rel(0.5), angle = 90)) 

print("Track gene -> DONE")

if(is.null(snps_track)){
  plots = c(plots,  annotation = genes)
} else {
  plots = c(plots, snps = snps_track,  annotation = genes)
}

pdf(paste0("done/",uniq_id,"_chipseq.pdf"))

selected_plots = create_plots_list(plots)

tracks(selected_plots, title = fig_title) +  xlim(current_range)
#tracks(c(plots[1:6], plots[19]), label.text.cex = 0.5)
#tracks(c(plots[7:12], plots[19]), label.text.cex = 0.5)
#tracks(plots[13:19], label.text.cex = 0.5)


dev.off()
save(means, file = paste0("done/means_",uniq_id,"_chipseq.rda"))
save(plots, file = paste0("done/plots_",uniq_id,"_chipseq_plots.rda"))

print(Sys.time() - t0)
