# SHINYSNP server.R

source("data_functions.R")

# EXECUTER 1 FOIS AU LANCEMENT DE LAPPLI
# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)
con = ""
logfile = tempfile(pattern = "log", tmpdir = "logs", fileext = "")
tempid = strsplit(x = logfile, split = "logs/log", fixed = TRUE)[[1]][2]
snpsfile = paste0("logs/snp", tempid)
hgsfile = paste0("logs/hg", tempid)
USER = "Audrey Lemacon"
print("27950000-28735000")

if(!file.exists(logfile)) {
  file.create(logfile)
  logcon=file(logfile,open="w")
  writeLines(c(USER,format(Sys.time(), "%a %b %d %X %Y")), logcon)
}

if(!file.exists(snpsfile)) {
  file.create(snpsfile)
  snpcon=file(snpsfile,open="a+b")
  writeLines(paste("snp_id","start","end","metadata",sep = "\t"), snpcon)
}

if(!file.exists(hgsfile)) {
  file.create(hgsfile)
  hgcon=file(hgsfile,open="a+b")
  writeLines(paste("start","end","method",sep = "\t"), hgcon)
}

print(snpsfile)
options(scipen=999) # virer l'annotation scientifique

shinyServer(function(input, output, session) {
  
  #### CHANGE POSITION VALUE ACCORDING TO THE SELECTED CHROM
  observe({
    selected_chr = input$chr
    selected_chr_min = subset(chroms, chr == selected_chr)$start
    selected_chr_max = subset(chroms, chr == selected_chr)$end
    
    updateNumericInput(session, "position_min", value =  selected_chr_min, 
                       min = selected_chr_min, max = selected_chr_max - 1)
    
    updateNumericInput(session, "position_max", value =  selected_chr_max, 
                       min = selected_chr_min + 1, max = selected_chr_max)    
  })
  
  observe ({
    if(!is.na(input$position_min)) {
      
      if (input$addsnp == 0)
        return()
      
      isolate ({
        
        if(is.na(input$snp_position_min)) {
          snpstart <- 0
        } else {
          snpstart <- input$snp_position_min
        }
        
        if(is.na(input$snp_position_max)) {
          snpend <- 0
        } else {
          snpend <- input$snp_position_max
        }
        
        selected_chr_min <- input$position_min
        selected_chr_max <- input$position_max
        
        error_msg <- c()
        iserror <- FALSE
        
        print(paste("start", snpstart))
        print(paste("end", snpend))
        
        if(snpend < snpstart){
          error_msg <- "SNP end must be greater than SNP start"
          iserror <- TRUE
          createAlert(session, "alert1", "exampleAlert1", title = "Warning",
                      content = error_msg, append = TRUE, dismiss = TRUE)
        } else {
          closeAlert(session, "exampleAlert1")
        }
        
        if(snpend == 0 || snpstart == 0){
          error_msg <- "SNP positions are mandatory"
          iserror <- TRUE
          createAlert(session, "alert2", "exampleAlert2", title = "Warning",
                      content = error_msg, append = TRUE, dismiss = TRUE)
        } else {
          closeAlert(session, "exampleAlert2")
        }
        
        if(snpstart >= selected_chr_min && snpstart <= selected_chr_max &&
             snpend >= selected_chr_min && snpend <= selected_chr_max) {
          closeAlert(session, "exampleAlert3")          
        } else {
          error_msg <- "The variant must be inclued in the selected region"
          iserror <- TRUE
          createAlert(session, "alert3", "exampleAlert3", title = "Warning",
                      content = error_msg, append = TRUE, dismiss = TRUE)
        }
        
        if(!iserror) {
          writeLines(text = paste(input$snp_label,snpstart,snpend,"n",sep="\t"),
                     con = snpcon)
          output$addsnp_msg <- renderText({
            return(paste0("Adding SNP ",input$snp_label," position : [", snpstart,
                          "-",snpend, "]"))})
        }
      })
    }
  })
  
  output$highlight <- renderUI({
    list(
      div(style="form-inline",
          tags$label("Higlight zone", `for` = "highlight-zone"), 
          div(style="input-daterange input-group",
              tags$input(id = "hgstart", type = "number", value = "start",class="input-small"),
              tags$span("to",style = "input-group-addon"),
              tags$input(id = "hgend", type = "number", value = "end",class="input-small")
          )
      ),
      br(),
      div(style="form-group shiny-input-container",
          tags$label("Higlight method", `for` = "highlight-method"), 
          checkboxGroupInput("hgmethod", label = NULL, 
                             choices = list("alpha" = "alpha", "color" = "color"),
                             inline = TRUE, selected = "alpha")
      ),
      actionButton(inputId = "addhg", label = "Add new Highlight Zone"),
      br(),
      br(),
      span(textOutput("addhg_msg"), style = "color:green")
    )
  })
  
  #   outputOptions(output, "addhg_msg", suspendWhenHidden=FALSE)
  
  output$addhg_msg <- renderText({
    if(!is.na(input$position_min)) {
      
      if (input$addhg == 0)
        return()
      
      isolate ({
        if(is.na(input$hgstart)) {
          hgstart <- 0
        } else {
          hgstart <- input$hgstart
        }
        
        if(is.na(input$hgend)) {
          hgend <- 0
        } else {
          hgend <- input$hgend
        }
        
        selected_chr_min <- input$position_min
        selected_chr_max <- input$position_max
        
        error_msg <- c()
        iserror <- FALSE
        
        if(hgend < hgstart){
          error_msg <- "HG Zone end must be greater than HG Zone start"
          iserror <- TRUE
          createAlert(session, "alert4", "exampleAlert4", title = "Warning",
                      content = error_msg, append = TRUE, dismiss = TRUE)
        } else {
          closeAlert(session, "exampleAlert4")
        }
        
        if(hgend == 0 || hgstart == 0){
          error_msg <- "HG Zone positions are mandatory"
          iserror <- TRUE
          createAlert(session, "alert5", "exampleAlert5", title = "Warning",
                      content = error_msg, append = TRUE, dismiss = TRUE)
        } else {
          closeAlert(session, "exampleAlert5")
        }
        
        if(hgstart >= selected_chr_min && hgstart <= selected_chr_max &&
             hgend >= selected_chr_min && hgend <= selected_chr_max) {
          closeAlert(session, "exampleAlert6")          
        } else {
          error_msg <- "The HG Zone must be inclued in the selected region"
          iserror <- TRUE
          createAlert(session, "alert6", "exampleAlert6", title = "Warning",
                      content = error_msg, append = TRUE, dismiss = TRUE)
        }
        
        if(!iserror) {
          writeLines(text = paste(hgstart,hgend,paste(input$hgmethod, collapse = ","),sep="\t"),
                     con = hgcon)
          return(paste0("Adding Highlight zone at position [", hgstart,
                        "-",hgend, "] with method : ",paste(input$hgmethod, collapse = ",")))
        }
      })
    }
  })
  
  output$plot1 <- renderPlot({
    if (input$draw == 0)
      return()
    
    isolate ({
      my.data <- load3DData(current_chr = input$chr)
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      my.ranges <- get3DDataOverview(my.data = my.data, current_range = current_range)
      
      archs_tracks <- drawArchs(ranges_list = my.ranges,
                                current_range = current_range,
                                highlight_range_list = NULL)
      
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      ### READ SNP FILE THEN CONVERT IN DATAFRAME
      lines = unlist(strsplit(readLines(con = snpcon), split = "\t", 
                              fixed = TRUE))
      snpsdf = data.frame(matrix(ncol = 4,data = lines[5:length(lines)], 
                                 byrow = TRUE))
      colnames(snpsdf) <- lines[1:4]
      
      snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
      
      t = c(archs_tracks, snp_track, annot_track)
      
      warning("DONE calculate the tracks",call. = FALSE)
      
      tracks(t) + xlim(current_range)
    })
    
  })
  
  outputOptions(output, "plot1", suspendWhenHidden=FALSE)
  
  output$plot2 <- renderPlot({
    if (input$draw == 0)
      return()
    
    isolate ({
      my.data <- load1DData(current_chr = input$chr)
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
#       my.ranges <- get3DDataOverview(my.data = my.data, current_range = current_range)
#       
#       archs_tracks <- drawArchs(ranges_list = my.ranges$arch, 
#                                 current_range = current_range,
#                                 highlight_range_list = NULL)
#       
#       annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
#       
#       ### READ SNP FILE THEN CONVERT IN DATAFRAME
#       lines = unlist(strsplit(readLines(con = snpcon), split = "\t", 
#                               fixed = TRUE))
#       snpsdf = data.frame(matrix(ncol = 4,data = lines[5:length(lines)], 
#                                  byrow = TRUE))
#       colnames(snpsdf) <- lines[1:4]
#       
#       snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
#       
#       t = c(archs_tracks, snp_track, annot_track)
#       
#       warning("DONE calculate the tracks",call. = FALSE)
#       
#       tracks(t) + xlim(current_range)
    })
    
  })
  
  outputOptions(output, "plot2", suspendWhenHidden=FALSE)
  
  
  
  
  observe({
    roots = c(wd='.')
    shinyFileChoose(input, 'loadsnp', session=session, roots=roots,
                    filetypes=c('','csv','txt'))
    pfile = parseFilePaths(roots, input$loadsnp)
    updateTextInput(session, "path",  value = pfile$datapath)
  })
  
  observe({
    if(file.exists(input$path)) {
      # get data from imported snp file
      imported_snp = read.table(input$path, header=TRUE, 
                                stringsAsFactors=FALSE, quote = "\"", sep="\t")
      for (i in (1:nrow(imported_snp))) {
        new_snp = imported_snp[i,]
        writeLines(text = paste(new_snp$snp_id,new_snp$start,new_snp$end,
                                new_snp$metadata,sep="\t"), con = snpcon)
      }
      
      output$addsnp_msg <- renderText({
        paste0("Importing following snps : ", paste(imported_snp$snp_id, 
                                                    collapse = ","))
      })
    }
    
  })
  
  observe({
    if(input$endAnalysis > 0) {
      updateButton(session, "endAnalysis", disabled = TRUE)
      if(!is.null(snpcon) && isOpen(snpcon)){close(snpcon)}
      if(!is.null(hgcon) && isOpen(hgcon)){close(hgcon)}
      if(!is.null(logcon) && isOpen(logcon)){
        writeLines(c(format(Sys.time(), "%a %b %d %X %Y")), logcon)
        close(logcon)
      }
      
      stopApp()
      #if(file.exists(snpsfile)) {file.remove(snpsfile)}
    }
  })
  
  
  # subset(sgp_k562_lane13, (sgp_k562_lane13$V2 >= current_start & sgp_k562_lane13$V2 <= current_stop) | (sgp_k562_lane13$V3 >= current_start & sgp_k562_lane13$V3 <= current_stop))
})

