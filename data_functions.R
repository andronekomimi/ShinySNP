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
  
  ### impet
  k562_impet = tryCatch({h5read(file = "data/data.h5",name = paste0("externe/K562/",current_chr,"/IM-PET/df"))}, error = function(e) {return(data.frame())} )
  mcf7_impet = tryCatch({h5read(file = "data/data.h5",name = paste0("externe/MCF7/",current_chr,"/IM-PET/df"))}, error = function(e) {return(data.frame())} )
  hmec_impet = tryCatch({h5read(file = "data/data.h5",name = paste0("externe/HMEC/",current_chr,"/IM-PET/df"))}, error = function(e) {return(data.frame())} )
  
  my.data = list(
    lncrna = lncrna, 
    lncrna_expr = lncrna_expr, 
    enhancers = enhancers,
    k562_super_enhancers = k562_super_enhancers, 
    mcf7_super_enhancers = mcf7_super_enhancers, 
    hmec_super_enhancers = hmec_super_enhancers,
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
  
  g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, 
                      color = rep(label,nrow(S)),
                      alpha = rep(0.5,nrow(S)))
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
                      color = rep(label,nrow(S)),
                      alpha = rep(0.5,nrow(S)))
  cat(paste0("-> Find ",nrow(S)," interaction(s) for ", label, "\n"))
  invisible(g_ranges)
}


create_gr_from_df <- function(current_range, my.df, label) {
  if(nrow(my.df) > 0){
    tmp = subsetData(my.df, current_range)
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
    tmp = subsetData(my.df, current_range)
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

get1DDataOverview <- function(my.data, current_range) {
  
  my.ranges = list(
    enhancers = create_gr_from_df(my.df = my.data$enhancers,current_range = current_range, label = "ENHANCERS"), 
    k562_super_enhancers = create_gr_from_df(my.df = my.data$k562_super_enhancers,current_range = current_range, label = "K562 SUPER ENHANCERS"),
    mcf7_super_enhancers = create_gr_from_df(my.df = my.data$mcf7_super_enhancers,current_range = current_range, label = "MCF7 SUPER ENHANCERS"), 
    hmec_super_enhancers = create_gr_from_df(my.df = my.data$hmec_super_enhancers,current_range = current_range, label = "HMEC SUPER ENHANCERS"), 
    k562_impet = create_gr_from_df(my.df = my.data$k562_impet,current_range = current_range, label = "K562 IM-PET"),
    mcf7_impet = create_gr_from_df(my.df = my.data$mcf7_impet,current_range = current_range, label = "MCF7 IM-PET"), 
    hmec_impet = create_gr_from_df(my.df = my.data$hmec_impet,current_range = current_range, label = "HMEC IM-PET"))
  
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
  
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    range = ranges_list[[i]]
    
    if(length(range) > 0) {
      track_title = gsub(x = range[1]$color, pattern = " ", replacement = "\n")
      range_track =  ggplot(data = range) + 
        geom_segment(size = 1) + ylab(track_title) +
        theme_bw() + xlim(current_range) + guides(color= TRUE) +
        theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
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
  
  snp_ids <- as.character(snps_df$snp_id)
  
  snps_ranges <- IRanges(as.numeric(as.character(snps_df$start)), as.numeric(as.character(snps_df$end)))
  snps <- GRanges(seqnames = as.character(seqnames(current_range)), ranges = snps_ranges, imp = snps_df$metadata)
  snps$name <- snp_ids
  
  snps_track <- autoplot(snps, aes(color=imp)) +
    geom_text(aes(x = start, y = 1, label=name, angle = 90, vjust=-1), size = 1, color = "blue") +
    theme_bw() + xlim(current_range) + ylab(label) + guides(color= FALSE) + 
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
  
  snps_track
}

drawLNCRNAFigures <- function(my.df, current_range) {
  
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
    
    
    my.ranges = GRanges(seqnames = as.character(seqnames(current_range)), 
                        ranges = IRanges(left,right), id = id, 
                        color = rep(x = "LNCRNA", times = length(left)))
    
    list(
      lncrna_track = drawSegment(ranges_list = list(my.ranges), 
                                 current_range = current_range),
      lncrna_hist = (ggplot(lncrna_df, aes(x = cell, y = fpkm)) + 
                       geom_bar(stat = "identity") + facet_grid(. ~ ts))
    )
  } else {
    print("No TsInfos in the requested area")
  }
  
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