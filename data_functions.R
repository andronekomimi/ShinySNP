suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggbio))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(biovizBase))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(rhdf5))

load3DData <- function(current_chr, my.dataset) {
  
  ### chiapet
  #pol2
  k562_chiapet = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/K562/",current_chr,"/ChIA-PET/df"))}, error = function(e) {return(data.frame())} )
  mcf7_chiapet = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/MCF7/",current_chr,"/ChIA-PET/df"))}, error = function(e) {return(data.frame())} )
  hmec_chiapet = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/HMEC/",current_chr,"/ChIA-PET/df"))}, error = function(e) {return(data.frame())} )
  
  #ctcf
  k562_chiapet_ctcf = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/K562/",current_chr,"/ChIA-PET/CTCF/df"))}, error = function(e) {return(data.frame())} )
  mcf7_chiapet_ctcf = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/MCF7/",current_chr,"/ChIA-PET/CTCF/df"))}, error = function(e) {return(data.frame())} )
  
  #er
  mcf7_chiapet_er = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/MCF7/",current_chr,"/ChIA-PET/ER/df"))}, error = function(e) {return(data.frame())} )
  
  
  ### hic
  k562_hic = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/K562/",current_chr,"/Hi-C/df"))}, error = function(e) {return(data.frame())} )
  mcf7_hic = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/MCF7/",current_chr,"/Hi-C/df"))}, error = function(e) {return(data.frame())} )
  hmec_hic = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/HMEC/",current_chr,"/Hi-C/df"))}, error = function(e) {return(data.frame())} )
  
  ### 3C
  k562_3c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/K562/",current_chr,"/3C/df"))}, error = function(e) {return(data.frame())} )
  mcf7_3c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/MCF7/",current_chr,"/3C/df"))}, error = function(e) {return(data.frame())} )
  hmec_3c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/HMEC/",current_chr,"/3C/df"))}, error = function(e) {return(data.frame())} )
  
  ### 4C
  k562_4c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/K562/",current_chr,"/4C/df"))}, error = function(e) {return(data.frame())} )
  mcf7_4c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/MCF7/",current_chr,"/4C/df"))}, error = function(e) {return(data.frame())} )
  hmec_4c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/HMEC/",current_chr,"/4C/df"))}, error = function(e) {return(data.frame())} )
  
  ### 5C 
  k562_5c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/K562/",current_chr,"/5C/df"))}, error = function(e) {return(data.frame())} )
  mcf7_5c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/MCF7/",current_chr,"/5C/df"))}, error = function(e) {return(data.frame())} )
  hmec_5c = tryCatch({h5read(file = "data/data.h5",name = paste0(my.dataset,"/HMEC/",current_chr,"/5C/df"))}, error = function(e) {return(data.frame())} )
  
  my.data = list(
    mcf7_chiapet = mcf7_chiapet, 
    hmec_chiapet = hmec_chiapet, 
    k562_chiapet = k562_chiapet,
    k562_chiapet_ctcf = k562_chiapet_ctcf,
    mcf7_chiapet_ctcf = mcf7_chiapet_ctcf,
    mcf7_chiapet_er = mcf7_chiapet_er,
    mcf7_hic = mcf7_hic, 
    hmec_hic = hmec_hic, 
    k562_hic = k562_hic,
    mcf7_3c = mcf7_3c, 
    hmec_3c = hmec_3c, 
    k562_3c = k562_3c,
    mcf7_4c = mcf7_4c, 
    hmec_4c = hmec_4c, 
    k562_4c = k562_4c,
    mcf7_5c = mcf7_5c, 
    hmec_5c = hmec_5c, 
    k562_5c = k562_5c
  )
  
  invisible(my.data)
  
}

load1DData <- function(current_chr) {
  ### lncrna 
  lncrna = tryCatch({h5read(file = "data/data.h5",name = paste0("common/",current_chr,"/LNCRNA/df"))}, error = function(e) {return(data.frame())} )
  
  ### lncrna_expr 
  lncrna_expr = tryCatch({h5read(file = "data/data.h5",name = paste0("common/misc/LNCRNA-EXP/df"))}, error = function(e) {return(data.frame())} )
  
  ### enhancers
  enhancers = tryCatch({h5read(file = "data/data.h5",name = paste0("common/",current_chr,"/ENHANCER/df"))}, error = function(e) {return(data.frame())} )
  
  ### super-enhancers
  k562_super_enhancers = tryCatch({h5read(file = "data/data.h5",name = paste0("common/K562/",current_chr,"/SUPER-ENHANCER/df"))}, error = function(e) {return(data.frame())} )
  mcf7_super_enhancers = tryCatch({h5read(file = "data/data.h5",name = paste0("common/MCF7/",current_chr,"/SUPER-ENHANCER/df"))}, error = function(e) {return(data.frame())} )
  hmec_super_enhancers = tryCatch({h5read(file = "data/data.h5",name = paste0("common/HMEC/",current_chr,"/SUPER-ENHANCER/df"))}, error = function(e) {return(data.frame())} )
  
  my.data = list(
    lncrna = lncrna, 
    lncrna_expr = lncrna_expr, 
    enhancers = enhancers,
    k562_super_enhancers = k562_super_enhancers, 
    mcf7_super_enhancers = mcf7_super_enhancers, 
    hmec_super_enhancers = hmec_super_enhancers
  )
  
  invisible(my.data)
  
  
}

load2DData <- function(current_chr) {
  
  ### impet
  k562_impet = tryCatch({h5read(file = "data/data.h5",name = paste0("externe/K562/",current_chr,"/IM-PET/df"))}, error = function(e) {return(data.frame())} )
  mcf7_impet = tryCatch({h5read(file = "data/data.h5",name = paste0("externe/MCF7/",current_chr,"/IM-PET/df"))}, error = function(e) {return(data.frame())} )
  hmec_impet = tryCatch({h5read(file = "data/data.h5",name = paste0("externe/HMEC/",current_chr,"/IM-PET/df"))}, error = function(e) {return(data.frame())} )
  
  my.data = list(
    k562_impet = k562_impet, 
    mcf7_impet = mcf7_impet, 
    hmec_impet = hmec_impet
  )
  
  invisible(my.data)
  
  
}

subsetData <- function(M, current_range){
  subset(M, (M$InteractorAChr == as.character(seqnames(current_range)) & 
               M$InteractorAStart >= start(current_range) & 
               M$InteractorAStart <= end(current_range)) | 
           (M$InteractorBEnd >= start(current_range) & 
              M$InteractorBEnd <= end(current_range)))
}

convert2GRange <- function(S, current_chr, label) {
  
  bool = (S$InteractorAStart < S$InteractorBEnd)
  
  left = c(S$InteractorAStart[bool], S$InteractorBEnd[!bool])
  right = c(S$InteractorBEnd[bool], S$InteractorAStart[!bool])
  
  ranges = IRanges(left,right)
  
  confidence = S$Confidence_Score1
  if(is.null(confidence)) {
    confidence = rep(0.5,nrow(S))
  }
  
  gene = S$Bgene
  if(is.null(gene)) {
    gene = rep("",nrow(S))
  }
  
  g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, 
                      label = rep(label,nrow(S)),
                      color = rep("black",nrow(S)),
                      gene = gene,
                      alpha = (1-confidence))
  cat(paste0("-> Find ",nrow(S)," interaction(s) for ", label, "\n"))
  invisible(g_ranges)
}


convert2GRange4IMPET <- function(S, current_chr, label) {
  is.gene = grepl(x = S$Bgene, pattern = "*.,ENSG.*")
  
  left = S$InteractorAStart[is.gene]
  right = S$InteractorAEnd[is.gene]
  confidence = S$Confidence_Score1[is.gene]
  gene = S$Bgene[is.gene]
  gene = gsub(x = gene, pattern = ",.*",replacement = "")
  dest = (S$InteractorBStart[is.gene] + S$InteractorBEnd[is.gene])/2
  
  df.temp = data.frame(left,right,gene,dest,confidence=(1-confidence))
  df.temp = unique(df.temp)
  
  ranges = IRanges(df.temp$left,df.temp$right)
  
  g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, 
                      label = rep(label,nrow(df.temp)),
                      color = rep("black",nrow(df.temp)),
                      gene = df.temp$gene,
                      dest = dest,
                      alpha = df.temp$confidence)
  
  cat(paste0("-> Find ",nrow(S)," interaction(s) for ", label, "\n"))
  
  invisible(g_ranges)
}

convert2GRange4LNCRNA <- function(S, current_chr, label) {
  
  bool = (S$InteractorAStart < S$InteractorBEnd)
  
  left = c(S$InteractorAStart[bool], S$InteractorBEnd[!bool])
  right = c(S$InteractorBEnd[bool], S$InteractorAStart[!bool])
  
  ranges = IRanges(left,right)
  
  g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, 
                      id = S$transcript_id,
                      color = rep("black",nrow(S)),
                      label = rep(label,nrow(S)),
                      alpha = rep(0.5,nrow(S)))
  cat(paste0("-> Find ",nrow(S)," interaction(s) for ", label, "\n"))
  invisible(g_ranges)
}



create_gr_from_df <- function(current_range, my.df, label) {
  if(nrow(my.df) > 0){
    tmp = unique(subsetData(my.df, current_range))
    if(nrow(tmp) > 0) {
      
      convert2GRange(tmp, current_chr = as.character(seqnames(current_range)), 
                     label = label)
    } else {
      invisible(GRanges())
    }
  } else  {
    invisible(GRanges())
  }
  
}

create_gr_from_df_4LNCRNA <- function(current_range, my.df, label) {
  if(nrow(my.df) > 0){
    tmp = unique(subsetData(my.df, current_range))
    if(nrow(tmp) > 0) {
      convert2GRange4LNCRNA(tmp, current_chr = as.character(seqnames(current_range)), 
                            label = label)
    } else {
      invisible(GRanges())
    }
  } else  {
    invisible(GRanges())
  }
  
}

create_gr_from_df_4IMPET <- function(current_range, my.df, label) {
  if(nrow(my.df) > 0){
    tmp = unique(subsetData(my.df, current_range))
    if(nrow(tmp) > 0) {
      has.gene = sum(grepl(x = tmp$Bgene, pattern = "*.,ENSG.*")) > 0
      if (has.gene){
        convert2GRange4IMPET(tmp, current_chr = as.character(seqnames(current_range)), 
                             label = label)
      } else {
        invisible(GRanges())
      }
    } else {
      invisible(GRanges())
    }
  } else  {
    invisible(GRanges())
  }
  
}


get1DDataOverview <- function(my.data, current_range) {
  
  my.ranges = list(
    enhancers = create_gr_from_df(my.df = my.data$enhancers,current_range = current_range, label = "ENHANCERS"), 
    k562_super_enhancers = create_gr_from_df(my.df = my.data$k562_super_enhancers,current_range = current_range, label = "K562 SUPER ENHANCERS"),
    mcf7_super_enhancers = create_gr_from_df(my.df = my.data$mcf7_super_enhancers,current_range = current_range, label = "MCF7 SUPER ENHANCERS"), 
    hmec_super_enhancers = create_gr_from_df(my.df = my.data$hmec_super_enhancers,current_range = current_range, label = "HMEC SUPER ENHANCERS"))
  
  invisible(my.ranges)
}

get2DDataOverview <- function(my.data, current_range) {
  
  my.ranges = list(
    k562_impet = create_gr_from_df_4IMPET(my.df = my.data$k562_impet,current_range = current_range, label = "K562 IM-PET"),
    mcf7_impet = create_gr_from_df_4IMPET(my.df = my.data$mcf7_impet,current_range = current_range, label = "MCF7 IM-PET"), 
    hmec_impet = create_gr_from_df_4IMPET(my.df = my.data$hmec_impet,current_range = current_range, label = "HMEC IM-PET"))
  
  invisible(my.ranges)
}

get3DDataOverview <- function(my.data, current_range) {
  
  my.ranges = list(
    mcf7_chiapet = create_gr_from_df(my.df = my.data$mcf7_chiapet,current_range = current_range, label = "MCF7 ChIA-PET"), 
    hmec_chiapet = create_gr_from_df(my.df = my.data$hmec_chiapet,current_range = current_range, label = "HMEC ChIA-PET"),
    k562_chiapet = create_gr_from_df(my.df = my.data$k562_chiapet,current_range = current_range, label = "K562 ChIA-PET"),
    k562_chiapet_ctcf = create_gr_from_df(my.df = my.data$k562_chiapet_ctcf,current_range = current_range, label = "K562 ChIA-PET CTCF"),
    mcf7_chiapet_ctcf = create_gr_from_df(my.df = my.data$mcf7_chiapet_ctcf,current_range = current_range, label = "MCF7 ChIA-PET CTCF"), 
    mcf7_chiapet_er = create_gr_from_df(my.df = my.data$mcf7_chiapet_er,current_range = current_range, label = "MCF7 ChIA-PET ER"), 
    mcf7_hic = create_gr_from_df(my.df = my.data$mcf7_hic,current_range = current_range, label = "MCF7 Hi-C"), 
    hmec_hic = create_gr_from_df(my.df = my.data$hmec_hic,current_range = current_range, label = "HMEC Hi-C"), 
    k562_hic = create_gr_from_df(my.df = my.data$k562_hic,current_range = current_range, label = "K562 Hi-C"),
    mcf7_3c = create_gr_from_df(my.df = my.data$mcf7_3c,current_range = current_range, label = "MCF7 3C"), 
    hmec_3c = create_gr_from_df(my.df = my.data$hmec_3c,current_range = current_range, label = "HMEC 3C"), 
    k562_3c = create_gr_from_df(my.df = my.data$k562_3c,current_range = current_range, label = "K562 3C"),
    mcf7_4c = create_gr_from_df(my.df = my.data$mcf7_4c,current_range = current_range, label = "MCF7 4C"), 
    hmec_4c = create_gr_from_df(my.df = my.data$hmec_4c,current_range = current_range, label = "HMEC 4C"), 
    k562_4c = create_gr_from_df(my.df = my.data$k562_4c,current_range = current_range, label = "K562 4C"),
    mcf7_5c = create_gr_from_df(my.df = my.data$mcf7_5c,current_range = current_range, label = "MCF7 5C"), 
    hmec_5c = create_gr_from_df(my.df = my.data$hmec_5c,current_range = current_range, label = "HMEC 5C"), 
    k562_5c = create_gr_from_df(my.df = my.data$k562_5c,current_range = current_range, label = "K562 5C"))
  
  invisible(my.ranges)
}

drawArchs <- function(ranges_list, highlight_ranges, current_range) {
  #highlight_range_list : start, stop, highlight method
  
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    my.range = ranges_list[[i]]
    if(length(my.range) > 0) {
      if(!is.null(highlight_ranges)) {
        my.range = make_emphasis(my.range, highlight_ranges)
      }
      
      track_title = gsub(x = my.range[1]$label, pattern = " ", replacement = "\n")
      range_track = get_chiapet_arch(my.range,track_title, current_range)
      my.tracks = c(my.tracks, range_track)
    }   
  }
  
  invisible(my.tracks)
}


getOverlaps <- function(ranges_list, snps_df, current_range) {
  
  overlaps = ""
  
  ### SNP Genomic Ranges
  ids <- as.character(snps_df$id)       
  snps_ranges <- IRanges(as.numeric(as.character(snps_df$start)), as.numeric(as.character(snps_df$end)))
  snps <- GRanges(seqnames = as.character(seqnames(current_range)), ranges = snps_ranges, imp = snps_df$metadata)
  snps$name <- ids
  
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    
    my.range = ranges_list[[i]]
    
    if(length(my.range) > 0) {
      range_name <-  my.range[1]$label
      overlaps <- paste0(overlaps, "Overlaps between variants and ",range_name,":\n")
      
      hits <- findOverlaps(query = my.range, subject = snps)
      
      for(i in seq_along(hits)) {
        hit = hits[i]
        query_idx = queryHits(hit)
        subject_idx = subjectHits(hit)
        overlap_snp <- snps[subject_idx]
        overlap_region <- my.range[query_idx]
        overlaps <- paste0(overlaps, "\t\t- ",overlap_snp$name, " (",start(overlap_snp),") <---> [",start(overlap_region),"-",end(overlap_region), "]\n")
      }
    }   
  }
  
  invisible(overlaps)
}

drawSegment <- function(ranges_list, highlight_ranges, current_range) {
  
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    my.range = ranges_list[[i]]
    
    if(length(my.range) > 0) {
      if(!is.null(highlight_ranges)) {
        my.range = make_emphasis(my.range, highlight_ranges)      
      }
      
      track_title = gsub(x = my.range[1]$label, pattern = " ", replacement = "\n")
      my.colors = unique(my.range$color)
      names(my.colors) = my.colors
      
      range_track =  ggbio::autoplot(my.range, aes(fill=color, alpha=alpha)) +
        scale_fill_manual(values = my.colors) +
        theme_bw() + xlim(current_range) + guides(fill= FALSE,alpha=FALSE) + 
        ylab(track_title) + theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
      
      if(!is.null(my.range$id)) {
        range_track = range_track + geom_text(aes(x = (start + ((end - start)/2)),
                                                  y = 1, label=id, angle=45), size = 2, color = "blue")
      }
      
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
    tryCatch({
      gr_txdb$symbols <- AnnotationDbi::select(org.Hs.eg.db,
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
    }, error = function(err) {
      print(err)
    })
    
  }
  
  return(g_track)
  
}

drawSNP <- function(current_range, snps_df, label) {
  
  ids <- as.character(snps_df$id)
  
  snps_ranges <- IRanges(as.numeric(as.character(snps_df$start)), as.numeric(as.character(snps_df$end)))
  snps <- GRanges(seqnames = as.character(seqnames(current_range)), ranges = snps_ranges, imp = snps_df$metadata)
  snps$name <- ids
  
  snps_track <- autoplot(snps, aes(fill=imp)) +
    geom_text(aes(x = start, y = 1, label=name, angle = 90, vjust=-1), size = 1, color = "blue") +
    theme_bw() + xlim(current_range) + ylab(label) + guides(fill= FALSE) + 
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
  
  snps_track
}

drawLNCRNAFigures <- function(my.df, current_range, highlight_ranges) {
  
  lncrna = create_gr_from_df_4LNCRNA(my.df,current_range = current_range, label = "LNCRNA")
  
  ts = unique(lncrna$id)
  df = h5read(file = "data/data.h5", name = "/common/misc/LNCRNA-EXP/df")
  ts_infos = unique(subset(df, df$transcript_id %in% ts))
  
  if(nrow(ts_infos) > 0) {
    mcf7_data = c()
    k562_data = c()
    hmec_data = ts_infos$encode_breast_hmec_cshl_rep1
    
    for(i in seq(1:nrow(ts_infos))) {
      mcf7_data = c(mcf7_data, mean(c(ts_infos$encode_breast_mcf7_caltech_rep1[i],
                                      ts_infos$encode_breast_mcf7_caltech_rep2[i],
                                      ts_infos$encode_breast_mcf7_caltech_rep3[i], 
                                      ts_infos$encode_breast_mcf7_cshl_rep1[i],
                                      ts_infos$encode_breast_mcf7_cshl_rep2[i])))
      
      k562_data = c(k562_data, mean(c(ts_infos$encode_cml_k562_cshl_rep1[i],
                                      ts_infos$encode_cml_k562_cshl_rep2[i], 
                                      ts_infos$encode_cml_k562_caltech_rep1[i],
                                      ts_infos$encode_cml_k562_caltech_rep2[i])))
      
    }
    
    lncrna_df = data.frame(ts =  rep(times = 3, x = ts_infos$transcript_id),
                           fpkm = c(k562_data, hmec_data, mcf7_data),
                           cell = c(rep(times = length(ts_infos$transcript_id), x = c("K562")), 
                                    rep(times = length(ts_infos$transcript_id), x = c("HMEC")), 
                                    rep(times = length(ts_infos$transcript_id), x = c("MCF7"))))
    
    
    ######## LNCRNA TRACK
    
    left = c()
    right = c()
    id = ts_infos$transcript_id
    
    #     # ranges
    for(i in seq(1:nrow(ts_infos))) {
      current_id = id[i]
      left = c(left, as.numeric(subset(my.df,my.df$transcript_id ==  current_id, select = c(InteractorAStart)))[1])
      right = c(right, as.numeric(subset(my.df,my.df$transcript_id ==  current_id, select = c(InteractorBEnd)))[1])
      
    }
    
    
    my.dataframe = data.frame(InteractorAStart = left, 
                              InteractorAStart = right, 
                              id = id, stringsAsFactors = FALSE)
    
    my.range = GRanges(seqnames = as.character(seqnames(current_range)), 
                       ranges = IRanges(left,right), id = id,  
                       alpha = rep(0.5,times = length(left)),
                       color = rep(x = "black", times = length(left)),
                       label = rep(x = "LNCRNA", times = length(left)))
    
    list(
      lncrna_track = drawSegment(ranges_list = list(my.range), 
                                 current_range = current_range, 
                                 highlight_ranges = highlight_ranges),
      lncrna_hist = (ggplot(lncrna_df, aes(x = cell, y = fpkm)) + 
                       geom_bar(stat = "identity") + 
                       geom_text(aes(label=round(fpkm, digits = 3),vjust=-1.3), size = 2, color = "blue") +
                       facet_grid(. ~ ts))
    )
  } else {
    print("No TsInfos in the requested area")
    df = data.frame(x=c(1), 
                    y=c(1), 
                    name = c("No annotated LNCRNA in this area"))
    
    list(
      lncrna_track = NULL,
      lncrna_hist = ggplot(data=df, mapping=aes(x=x, y=y)) +
        geom_blank() + ylab("") + xlab("") + 
        geom_text(aes(x = x, y = y, label=name), size = 7, color = "red") +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
    )
  }
  
}

drawIMPET <- function(ranges_list, highlight_ranges, current_range){
  
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    my.range = ranges_list[[i]]
    
    if(length(my.range) > 0) {
      #       if(!is.null(highlight_ranges)) {
      #         my.range = make_emphasis(my.range, highlight_ranges)      
      #       }
      
      df_impet <- data.frame(x = (start(my.range)), xend = (my.range$dest), gene = my.range$gene)
      
      track_title = gsub(x = my.range[1]$label, pattern = " ", replacement = "\n")
      my.colors = unique(my.range$color)
      names(my.colors) = my.colors
      
      g = autoplot(my.range, aes(color=gene, fill = gene, group = gene, alpha = alpha), group.selfish = FALSE) +
        geom_arch(df_impet, aes(x = x, xend=xend,color = gene, group = gene))
      
      range_track = g + theme_bw() + xlim(current_range) + 
        guides(alpha=FALSE, col = guide_legend(ncol =  2, keywidth = 0.2, keyheight = 0.4), 
               fill = guide_legend(ncol =  2, keywidth = 0.2, keyheight = 0.4)) + 
        ylab(track_title) + theme(axis.title.y = element_text(size = rel(0.5), angle = 90),
                                  legend.text =  element_text(size = rel(0.4)),
                                  legend.title = element_text(size = rel(0.4)))
      
      if(!is.null(my.range$id)) {
        range_track = range_track + geom_text(aes(x = (start + ((end - start)/2)),
                                                  y = 1, label=id, angle=45), size = 2, color = "blue")
      }
      
      my.tracks = c(my.tracks, range_track)
    }   
  }
  
  invisible(my.tracks)
}



mergeRanges <- function(ranges, label) {
  can_run_6(label)
  ranges$label = rep(x = label, times = length(ranges))
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
  my.colors = unique(rep$color)
  names(my.colors) = my.colors
  
  g = ggplot(rep) +
    geom_arch() +
    theme_bw() +
    aes(color=color, alpha = alpha) +
    #   aes(color=color) +
    #scale_colour_manual(values = c("gray","promoter"="black")) +
    scale_colour_manual(values = my.colors) +
    xlim(current_range) +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
    ylab(track_title) + guides(alpha=FALSE, color=FALSE) +
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
  
  g
  
}

make_emphasis = function(my.ranges, hg.ranges) {
  hits = findOverlaps(query = my.ranges, subject = hg.ranges)
  
  for(i in seq_along(hits)) {
    hit = hits[i]
    query_idx = queryHits(hit)
    subject_idx = subjectHits(hit)
    
    #my.ranges[query_idx]$alpha = hg.ranges[subject_idx]$alpha
    my.ranges[query_idx]$color = hg.ranges[subject_idx]$color
  }
  
  my.ranges
}

#######################         RNASEQ/CHIPSEQ               ########################
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

get_views <- function(coverage, gr) { 
  chr <- intersect(names(coverage), as.character(seqlevels(gr)))
  Views(coverage[chr], as(gr, "RangesList")[chr])
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

drawRNASEQ <- function(file_list, highlight_file, current_range) {
  views_dir <- read.csv(file = "/etc/shiny-apps/ShinySNP.conf",header = TRUE)$VIEWS_DIR
  nb.cores <- as.numeric(read.csv(file = "/etc/shiny-apps/ShinySNP.conf",header = TRUE)$CORES)
  
  current_start <- start(current_range)
  current_stop <- end(current_range)
  current_chr <- as.character(seqnames(current_range))
  
  ret <- parse_file_list(file_list)
  bam_files_list <- ret[[1]]
  controls <- ret[[2]]
  bam_files <- unlist(bam_files_list)
  
  coverages <- list()
  for(n in names(bam_files_list)) {
    replicats = bam_files_list[[n]]
    print(replicats)
    means <- c()
    
    for(replicat in replicats) {
      load(paste0(views_dir,replicat,"_",current_chr,".Rda")) #views
      means <- colSums(rbind(means, as.vector(views$coverages[[current_chr]][[1]][current_start:current_stop])), na.rm=TRUE)
      remove(views)
    }
    
    means <- means / (length(replicats))
    coverages[[n]] <- means
    remove(means)
  }
  
  top_value <- max(unlist(mclapply(coverages, max, mc.cores = nb.cores)))
  trackname <- character()
  
  get_plot <- function(bam_file_name) {
    coverage <- coverages[[bam_file_name]]
    position <- seq(start(current_range), end(current_range))
    trackname <- bam_file_name
    data <- data.frame(position = position, coverage = coverage)
    
    if(!is.null(highlight_file) && file.exists(highlight_file)) {
      hgs_df <- read.table(highlight_file, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      #d <- data.frame(x1=as.numeric(hgs_df$start), x2=as.numeric(hgs_df$end), y1=0, y2=max(data$coverage), col = hgs_df$color)
      if(nrow(hgs_df) > 0){
        d <- data.frame(x1=as.numeric(hgs_df$start), x2=as.numeric(hgs_df$end), y1=0, y2=top_value, col = hgs_df$color)
      } else {
        d <- data.frame(x1=c(0), x2=c(0), y1=0, y2=top_value, col=c("black"))
      }
    } else {
      d <- data.frame(x1=c(0), x2=c(0), y1=0, y2=top_value, col=c("black"))
    }
    
    my.colors = d$col
    names(my.colors) = my.colors
    
    g = ggplot() + geom_line(data = data, mapping =  aes(x = position, y = coverage)) + 
      geom_rect(data = d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.4, fill = my.colors) +
      ylim(0,top_value) + ylab(trackname) + theme_bw() + xlim(current_range) + 
      theme(axis.title.y = element_text(angle = 0))
    
    print(paste0("Track plot for ", trackname," -> DONE"))
    
    g
  }
  
  plots <- mclapply(names(coverages), get_plot, mc.cores = nb.cores)
  remove(coverages)
  invisible(plots)
}


drawCHIPSEQ <- function(file_list, highlight_file, current_range) {
  
  views_dir <- read.csv(file = "/etc/shiny-apps/ShinySNP.conf",header = TRUE)$VIEWS_DIR
  nb.cores <- as.numeric(read.csv(file = "/etc/shiny-apps/ShinySNP.conf",header = TRUE)$CORES)
  
  current_start <- start(current_range)
  current_stop <- end(current_range)
  current_chr <- as.character(seqnames(current_range))
  
  ret <- parse_file_list(file_list)
  bam_files_list <- ret[[1]]
  controls <- ret[[2]]
  bam_files <- unlist(bam_files_list)
  
  coverages <- list()
  for(n in names(bam_files_list)) {
    replicats = bam_files_list[[n]]
    print(replicats)
    means <- c()
    
    for(replicat in replicats) {
      load(paste0(views_dir,replicat,"_",current_chr,".Rda")) #views
      means <- colSums(rbind(means, as.vector(views$coverages[[current_chr]][[1]][current_start:current_stop])), na.rm=TRUE)
      remove(views)
    }
    
    means <- means / (length(replicats))
    coverages[[n]] <- means
    remove(means)
  }
  
  top_value <- max(unlist(mclapply(coverages, max, mc.cores = nb.cores)))
  
  apply_control <- function(experiment) {
    ctr =  controls[[experiment]]
    if (! is.na(ctr)) {
      new_mean = coverages[[experiment]] - coverages[[ctr]]
    }
    else
    {
      new_mean = coverages[[experiment]]
    }
    
    new_mean
  }
  
  
  new_coverages = mclapply(names(coverages), apply_control, mc.cores = nb.cores)
  names(new_coverages) = names(coverages)
  coverages = new_coverages
  remove(new_coverages)
  
  trackname <- character()
  
  get_plot <- function(bam_file_name) {
    coverage <- coverages[[bam_file_name]]
    position <- seq(start(current_range), end(current_range))
    trackname <- bam_file_name
    data <- data.frame(position = position, coverage = coverage)
    
    if(!is.null(highlight_file) && file.exists(highlight_file)) {
      hgs_df <- read.table(highlight_file, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      #d <- data.frame(x1=as.numeric(hgs_df$start), x2=as.numeric(hgs_df$end), y1=0, y2=max(data$coverage), col = hgs_df$color)
      if(nrow(hgs_df) > 0){
        d <- data.frame(x1=as.numeric(hgs_df$start), x2=as.numeric(hgs_df$end), y1=0, y2=top_value, col = hgs_df$color)
      } else {
        d <- data.frame(x1=c(0), x2=c(0), y1=0, y2=top_value, col=c("black"))
      }
    } else {
      d <- data.frame(x1=c(0), x2=c(0), y1=0, y2=top_value, col=c("black"))
    }
    
    my.colors = d$col
    names(my.colors) = my.colors
    
    g = ggplot() + geom_line(data = data, mapping =  aes(x = position, y = coverage)) + 
      geom_rect(data = d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.4, fill = my.colors) +
      ylim(0,top_value) + ylab(trackname) + theme_bw() + xlim(current_range) +
      theme(axis.title.y = element_text(angle = 0))
    
    print(paste0("Track plot for ", trackname," -> DONE"))
    
    g
  }
  
  plots <- mclapply(names(coverages), get_plot, mc.cores = nb.cores)
  names(plots) = names(coverages)
  selected_plots = create_plots_list(plots)
  
  remove(coverages)
  remove(plots)
  
  invisible(selected_plots)
}


