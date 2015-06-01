suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggbio))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(biovizBase))
suppressMessages(library(org.Hs.eg.db))


load3DData <- function(current_chr) {
  
  ### chiapet 
  sgp_k562_lane13 <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane13.Rda"))
  sgp_k562_lane24 <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane24.Rda")) 
  sgp_mcf7_lane11 <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane11.Rda")) 
  sgp_mcf7_lane23 <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane23.Rda"))
  sgp_k562_lane13_inter <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane13_inter.Rda"))
  sgp_k562_lane24_inter <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane24_inter.Rda")) 
  sgp_mcf7_lane11_inter <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane11_inter.Rda"))
  sgp_mcf7_lane23_inter <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane23_inter.Rda"))
  sfd_k562_rep1 <- readRDS(paste0("data/",current_chr,"_sfd_k562_rep1.Rda"))  
  sfd_k562_rep2 <- readRDS(paste0("data/",current_chr,"_sfd_k562_rep2.Rda"))
  sfd_mcf7_rep3 <- readRDS(paste0("data/",current_chr,"_sfd_mcf7_rep3.Rda"))
  sfd_mcf7_rep4 <- readRDS(paste0("data/",current_chr,"_sfd_mcf7_rep4.Rda"))
  
  ### hic 
  hmec_file_1 <- readRDS(paste0("data/",current_chr,"_hmec_file_1.Rda"))
  hmec_file_2 <- readRDS(paste0("data/",current_chr,"_hmec_file_2.Rda"))
  k562_file_1 <- readRDS(paste0("data/",current_chr,"_k562_file_1.Rda"))
  k562_file_2 <- readRDS(paste0("data/",current_chr,"_k562_file_2.Rda"))
  
  my.data = list(sgp_k562_lane13 = sgp_k562_lane13,
                 sgp_k562_lane24 = sgp_k562_lane24,
                 sgp_mcf7_lane11 = sgp_mcf7_lane11,
                 sgp_mcf7_lane23 = sgp_mcf7_lane23,
                 sgp_k562_lane13_inter = sgp_k562_lane13_inter,
                 sgp_k562_lane24_inter = sgp_k562_lane24_inter,
                 sgp_mcf7_lane11_inter = sgp_mcf7_lane11_inter,
                 sgp_mcf7_lane23_inter = sgp_mcf7_lane23_inter,
                 sfd_k562_rep1 = sfd_k562_rep1,
                 sfd_k562_rep2 = sfd_k562_rep2,
                 sfd_mcf7_rep3 = sfd_mcf7_rep3,
                 sfd_mcf7_rep4 = sfd_mcf7_rep4,
                 hmec_file_1 = hmec_file_1,
                 hmec_file_2 = hmec_file_1,
                 k562_file_1 = k562_file_1,
                 k562_file_2 = k562_file_2
  )
  
  invisible(my.data)
  
}

load1Data <- function(current_chr) {
  
  ### lncrna
  lncrna <- readRDS(paste0("data/",current_chr,"_mitrans.Rda"))
  lncrna_expr <- readRDS(paste0("data/lncrna_expr.Rda"))
  
  my.data = list(
    lncrna = lncrna,
    lncrna_expr = lncrna_expr
  )
  
  invisible(my.data)
  
}


subsetData <- function(my.df, start_idx, end_idx, zone_start, zone_stop){
  subset(my.df, (my.df[start_idx] >= zone_start & my.df[start_idx] <= zone_stop) | 
           (my.df[end_idx] >= zone_start & my.df[end_idx] <= zone_stop), 
         select = c(start_idx, end_idx))
  
}

convert2GRange <- function(my.df, current_chr, label) {
  left = c()
  right = c()
  color = c()
  alpha = c()  
  
  for (i in (1:nrow(my.df)))
  {
    start_value = as.numeric(my.df[i,1])
    end_value = as.numeric(my.df[i,2])
    
    left = c(left, start_value)
    right = c(right, end_value)
    color = c(color,label)
  }
  
  ranges = IRanges(left,right)
  
  g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, color = color, alpha = rep(0.5,length(left)))
  cat(paste0("-> Find ",length(left)," interaction(s) for ", label, "\n"))
  invisible(g_ranges)
}


create_gr_from_df <- function(current_range, my.df,start_idx,stop_idx, label) {
  tmp = subsetData(my.df,start_idx,stop_idx, start(current_range),end(current_range))
  if(nrow(tmp) > 0) {
    convert2GRange(tmp, current_chr = as.character(seqnames(current_range)), 
                   label = label)
  } else {
    invisible(GRanges())
  }
  
  
}

get3DDataOverview <- function(my.data, current_range) {
  
  my.ranges = list(
    sgp_k562_lane13 = create_gr_from_df(current_range,my.data$sgp_k562_lane13,2,6,"K562 ChIA-Pet lane 13"),
    sgp_k562_lane24 = create_gr_from_df(current_range,my.data$sgp_k562_lane24, 2,6,"K562 ChIA-Pet lane 24"),
    sgp_mcf7_lane11 = create_gr_from_df(current_range,my.data$sgp_mcf7_lane11,2,6,"MCF7 ChIA-Pet lane 11"),
    sgp_mcf7_lane23 = create_gr_from_df(current_range,my.data$sgp_mcf7_lane23,2,6,"MCF7 ChIA-Pet lane 23"),
    sfd_k562_rep1 = create_gr_from_df(current_range,my.data$sfd_k562_rep1,2,6,"K562 ChIA-Pet rep 1"),
    sfd_k562_rep2 = create_gr_from_df(current_range,my.data$sfd_k562_rep2,2,6,"K562 ChIA-Pet rep 2"),
    sfd_mcf7_rep3 = create_gr_from_df(current_range,my.data$sfd_mcf7_rep3,2,6,"MCF7 ChIA-Pet rep 3"),
    sfd_mcf7_rep4 = create_gr_from_df(current_range,my.data$sfd_mcf7_rep4,2,6,"MCF7 ChIA-Pet rep 4"),
    hmec_file_1 = create_gr_from_df(current_range,my.data$hmec_file_1,2,3,"HMEC HiC-arrowhead"),
    hmec_file_2 = create_gr_from_df(current_range,my.data$hmec_file_2,2,6,"HMEC HiC-hiccups"),
    k562_file_1 = create_gr_from_df(current_range,my.data$k562_file_1,2,3,"K562 HiC-arrowhead"),
    k562_file_2 = create_gr_from_df(current_range,my.data$k562_file_2,2,6,"K562 HiC-hiccups"))
  
  invisible(my.ranges)
}


drawArchs <- function(ranges_list = NULL, highlight_range_list = NULL, current_range) {
  #highlight_range_list : start, stop, highlight method
  
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    range = ranges_list[[i]]
    if(length(range) > 0) {
      for( i in seq_along(highlight_range_list)) {
        highlight_range = highlight_range_list[i]
        range = make_emphasis(range, highlight_range)      
      }
      track_title = gsub(x = range[1]$color, pattern = " ", replacement = "\n")
      range_track = get_chiapet_arch(range,track_title, current_range)
      my.tracks = c(my.tracks, range_track)
    }   
  }
  
  invisible(my.tracks)
}

drawSegment <- function(ranges_list = NULL, current_range) {
  
  if(is.null(ranges_list))
  {
    ranges_list = my.ranges$ranges$lncrna
  }
  
  print(ranges_list)
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    range = ranges_list[[i]]
    if(length(range) > 0) {
      track_title = gsub(x = range[1]$color, pattern = " ", replacement = "\n")
      range_track =  ggplot(data = range) + 
        geom_segment(mapping=aes(x=start, xend=end, color=exp), size = 10) + ylab("") +
        geom_text(aes(x = start, y = 1, label=tsID, vjust=3), size = 2, color = "blue") + 
        theme_bw() + xlim(current_range) + guides(color= TRUE)
      my.tracks = c(my.tracks, range_track)
    }   
  }
  
  invisible(my.tracks)
}

setStudyRange <- function(current_chr, current_start, current_stop) {
  current_range <- IRanges::IRanges(current_start, current_stop)
  current_range <- GenomicRanges::GRanges(seqnames = current_chr, ranges = current_range)
  current_range
}


setHighLight <- function(current_start, current_stop, method) {
  can_run_3(current_start, current_stop, method)
  
  # method : alpha, color, size
  special_range <- IRanges(current_start, current_stop)
  special_range <- GRanges(seqnames = current_chr, ranges = special_range, 
                           method = method)
  special_range
}


drawAnnotations <- function(label = "Annotations", current_range) {
  #can_run_4()
  
  df = data.frame(x=c(mean(c(start(current_range), end(current_range)))), 
                  y=c(1), 
                  name = c("No genomic element in this area"))
  
  err_label <- "No genomic element in this area"
  label_position = mean(c(start(current_range), end(current_range)))
  range <- IRanges(label_position, label_position)
  range <- GRanges(seqnames = as.character(seqlevels(current_range)), ranges = range)
  range$name = err_label
  
  g_track <- ggplot(range) +
    geom_text(aes(x = start, y = 1, label=name), size = 5, color = "red") +
    theme_bw() + xlim(current_range) + ylab("") + xlab("") + guides(color= FALSE) + 
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())
  
  gr_txdb <- NULL
  
  tryCatch({
    gr_txdb <- crunch(TxDb.Hsapiens.UCSC.hg19.knownGene, which = current_range)
  }, error = function(err) {
    print(err)
  })
  
  if(!is.null(gr_txdb)){
    colnames(values(gr_txdb))[4] <- "model"
    gr_txdb$symbols <- select(org.Hs.eg.db,
                              keys = as.character(gr_txdb$gene_id),
                              column = "SYMBOL",
                              keytype = "ENTREZID")$SYMBOL
    i <- which(gr_txdb$model == "gap")
    
    gr_txdb <- gr_txdb[-i]
    levels(gr_txdb) <- c("cds", "exon", "utr")
    grl_txdb <- split(gr_txdb, gr_txdb$symbols)
    
    if(length(grl_txdb) > 0){
      g_track <- autoplot(grl_txdb, aes(type = model)) + 
        theme_bw() + xlim(current_range) + ylab(label) +
        theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
    }
  }
  
  return(g_track)
  
}

drawSNP <- function(current_range, snps_df, label) {
  
  snp_ids <- snps_df$snp_id
  
  snps_ranges <- IRanges(as.numeric(as.character(snps_df$start)), as.numeric(as.character(snps_df$end)))
  snps <- GRanges(seqnames = as.character(seqnames(current_range)), ranges = snps_ranges, imp = snps_df$metadata)
  snps$name <- snp_ids
  
  snps_track <- autoplot(snps, aes(color=imp)) +
    geom_text(aes(x = start, y = 1, label=name, angle = 90, vjust=-1), size = 1, color = "blue") +
    theme_bw() + xlim(current_range) + ylab(label) + guides(color= FALSE) + 
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
  
  snps_track
}


mergeRanges <- function(ranges, label) {
  can_run_6(label)
  ranges$color = rep(x = label, times = length(ranges))
  ranges
}


organizeRanges <- function(track_index = NULL, range_index = NULL, labels = NULL) {
  
  can_run_5(track_index, range_index, labels)
  
  selected_range = c()
  
  if(!is.null(range_index)) {
    
    str = strsplit(x = range_index, split = ";", fixed = TRUE)
    
    for(j in seq_along(str)) {
      
      idx = as.integer(str[[j]])
      if(length(idx) > 1) { # merge
        temp = c()
        for(i in idx) {
          temp = c(temp, my.ranges$arch[[i]])
        }
        
        if(is.na(labels[j])) {
          can_run_6(labels[j])
        }
        else
        {
          r = mergeRanges(do.call("c",temp), label = labels[j])
        }
        
      }
      else # simple selection
      {
        if(!is.na(labels[j])) {
          r = mergeRanges(my.ranges$arch[[idx]], label = labels[j])
        }
        else
        {
          r = my.ranges$arch[[idx]]
        }
      }
      
      selected_range = c(selected_range, r)
    }
  }
  else # not range_index so track_index
  {
    str = strsplit(x = track_index, split = ";", fixed = TRUE)
    
    temp.ranges = c()
    
    for (my.range in my.ranges$arch) {
      if(length(my.range) > 0) {
        temp.ranges = c(temp.ranges, my.range)
      }
    }
    
    
    for(j in seq_along(str)) {
      
      idx = as.integer(str[[j]])
      
      if(length(idx) > 1) { # merge
        temp = c()
        for(i in idx) {
          temp = c(temp, temp.ranges[[i]])
        }
        
        if(is.na(labels[j])) {
          can_run_6(labels[j])
        }
        else
        {
          r = mergeRanges(do.call("c",temp), label = labels[j])
        }
      }
      else # simple selection
      {
        if(!is.na(labels[j])) {
          r = mergeRanges(temp.ranges[[idx]], label = labels[j])
        }
        else
        {
          r = temp.ranges[[idx]]
        }
      }
      
      selected_range = c(selected_range, r)
    }
  }
  
  selected_range
}

can_run <- function() {
  if(!exists("current_chr"))
  {
    stop(call. = FALSE, "Please define the chromosome to study before running 
         this function(ex : current_chr = \"chr12\")")
  }
}

can_run_2 <- function() {
  if(!exists("current_range"))
  {
    stop(call. = FALSE, "Please define the range of the chromosome to study 
         (ex : current_range = setStudyRange(27950000, 28735000))")
  }
  
  if(!exists("my.data"))
  {
    stop(call. = FALSE, "Please load data before running this function
         ex : my.data = loadChrData()")
  }
}

can_run_3 <- function(current_start, current_stop, method) {
  if(current_start < start(current_range) || current_stop > end(current_range))
  {
    stop(call. = FALSE, "The highlight zone must be included in the study range")
  }
  
  if(!(method %in% c("color","alpha")))
  {
    stop(call. = FALSE, "The highlight method must be \'color\' or \'alpha\'")
  }
}

can_run_4 <- function() {
  if(!exists("current_range"))
  {
    stop(call. = FALSE, "Please define the range of the chromosome to study 
         (ex : current_range = setStudyRange(27950000, 28735000))")
  }
}

can_run_5 <- function(track_index, range_index, labels) {
  
  if (is.null(track_index) && is.null(range_index)) {
    stop(call. = FALSE, "Please provide an index vector")
  }
  
  if (!is.null(track_index) && !is.null(range_index)) {
    stop(call. = FALSE, "Please provide only ONE index vector")
  }
  
  if (is.null(labels)) {
    stop(call. = FALSE, "Please provide a label vector. Use NA value for each 
         track if you don't want to rename them")
  }
  
  if (length(track_index) != length(labels) && 
        length(range_index) != length(labels)) {
    stop(call. = FALSE, "Please provide a label for EACH track. Use NA value 
        for each track if you don't want to rename them")
  }
}

can_run_6 <- function(label) {
  
  if (is.null(label) || is.na(label)) {
    stop(call. = FALSE, "Please provide a label to name the merged range")
  }
}

get_chiapet_arch <- function(rep, track_title, current_range) {
  
  g = ggplot(rep) +
    geom_arch() +
    theme_bw() +
    aes(color=color, alpha = alpha) +
    #   aes(color=color) +
    #   scale_colour_manual(values = c("gray","promoter"="black")) +
    xlim(current_range) +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
    ylab(track_title) + guides(alpha=FALSE, color=FALSE) +
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
  
  g
  
}