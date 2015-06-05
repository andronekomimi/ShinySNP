# SHINYSNP server.R

source("data_functions.R")
library(sendmailR)

# EXECUTER 1 FOIS AU LANCEMENT DE LAPPLI
# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)
logfile = tempfile(pattern = "log_", tmpdir = "logs", fileext = "")
tempid = strsplit(x = logfile, split = "logs/log_", fixed = TRUE)[[1]][2]
snpsfile = paste0("logs/snp_", tempid)
hgsfile = paste0("logs/hg_", tempid)
todornaseqfile = paste0("todo/rnaseq_", tempid)
todochipseqfile = paste0("todo/chipseq_", tempid)

USER = "Audrey Lemacon"
USERMAIL = "audrey.lemacon.1@ulaval.ca"
OPERATOR = "audrey.lemacon.1@ulaval.ca"

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

if(!file.exists(todornaseqfile)) {
  file.create(todornaseqfile)
  todornaseqcon=file(todornaseqfile,open="w")
  writeLines(paste0("#", USER), todornaseqcon)
  writeLines(paste0("#", tempid), todornaseqcon)
}

if(!file.exists(todochipseqfile)) {
  file.create(todochipseqfile)
  todochipseqcon=file(todochipseqfile,open="w")
  writeLines(paste0("#", USER), todornaseqcon)
  writeLines(paste0("#", tempid), todornaseqcon)
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
  
  observe({
    
    if(input$runEx == 0)
      return()
    
    isolate({
      updateNumericInput(session, "position_min", value = 27950000)
      updateNumericInput(session, "position_max", value = 28735000)
      updateNumericInput(session, "snp_position_min", value = 28155080)
      updateNumericInput(session, "snp_position_max", value = 28155080)
      
    })
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
          # si snpadd alors desactivation loadsnp
          updateButton(session, "loadsnp", disabled = TRUE)
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
    if (input$run == 0)
      return()
    
    isolate ({
      if(!is.null(snpcon) && isOpen(snpcon, "a+b")){
        close(con = snpcon)
        # effacer le fichier
      }
      
      if(!is.null(hgcon) && isOpen(hgcon, "a+b")){
        close(hgcon)
        # effacer le fichier
      }
      
      my.data <- load3DData(current_chr = input$chr, my.dataset = input$my.dataset)
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      my.ranges <- get3DDataOverview(my.data = my.data, current_range = current_range)
      
      archs_tracks <- drawArchs(ranges_list = my.ranges,
                                current_range = current_range,
                                highlight_range_list = NULL)
      
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      ### READ SNP FILE THEN CONVERT IN DATAFRAME
      snpsdf = read.table(snpsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      print(snpsdf)
      
      if(nrow(snpsdf) > 0 ){
        snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
        t = c(archs_tracks, snp_track, annot_track)
      } else {
        t = c(archs_tracks, annot_track)
      }
      
      warning("DONE calculate the tracks",call. = FALSE)
      tracks(t) + xlim(current_range)
    })
    
    
  })
  
  outputOptions(output, "plot1", suspendWhenHidden=FALSE)
  
  
  
  output$plot2 <- renderPlot({
    if (input$run == 0)
      return()
    
    isolate ({
      
      my.data <- load1DData(current_chr = input$chr)
      
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      my.ranges <- get1DDataOverview(my.data = my.data, current_range = current_range)
      
      segments_tracks <- drawSegment(ranges_list = my.ranges, 
                                current_range = current_range)
      
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      ### READ SNP FILE THEN CONVERT IN DATAFRAME
      snpsdf = read.table(snpsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      
      print(snpsdf)
      
      if(nrow(snpsdf) > 0 ){
        snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
        t = c(segments_tracks, snp_track, annot_track)
      } else {
        t = c(segments_tracks, annot_track)
      }
      
      warning("DONE calculate the tracks",call. = FALSE)
      
      tracks(t) + xlim(current_range)
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
    
    if (is.null(input$loadsnp))
      return()
    
    isolate ({
      if (file.exists(input$path)) {
        
        print("craquage")
        output$addsnp_msg <- renderText({
          return("PARRCE")
        })
        
        output$addsnp_err <- renderText({
          return("QUOI ???")
        })
        
        print("total")
        # desactive le bouton addsnp
        updateButton(session, "addsnp", disabled = TRUE)
        
        tryCatch ({
          # get data from imported snp file
          imported_snp = read.table(input$path, header=TRUE, 
                                    stringsAsFactors=FALSE, quote = "\"", sep="\t")
          if(is.null(snpcon)) {
            print('ok')
          }
          print(imported_snp)
          
          if((!"snp_id" %in% names(imported_snp)) || (!"start" %in% names(imported_snp)) ||
               (!"end" %in% names(imported_snp)) || (!"metadata" %in% names(imported_snp))) {
            stop(paste0("File format error : waiting for a tab delimited file containing ",
                        "at least the following 4 columns : snp_id, start, end, metadata"))
          } else {
            for (i in (1:nrow(imported_snp))) {
              new_snp = imported_snp[i,]
              writeLines(text = paste(new_snp$snp_id,new_snp$start,new_snp$end,
                                      new_snp$metadata,sep="\t"), con = snpcon)
            }
            
            print(paste0("Importing following snps : ", paste(imported_snp$snp_id, 
                                                              collapse = ",")))
            
            updateButton(session = session, label = "RHOOOOO",inputId = "loadsnp")
            
            output$addsnp_msg <- renderText({
              return(paste0("Importing following snps : ", paste(imported_snp$snp_id, 
                                                                 collapse = ",")))})
            
          }
        },
        error = function(e) {
          print(e)
          output$addsnp_err <- renderText({
            return(paste0("Error while trying to import snps file : ", e))
          })
        }
        )
        
      }
    })
  })
  
  observe({
    
    if(input$run == 0)
      return()
    
    isolate({
      updateButton(session, "run", disabled = TRUE)
      
    })
  })
  
  observe ({
    
    if(input$addRNASeq == 0) {
      return()
    }
    
    isolate ({
      # working files
      all_cell_merge_file = "all_cells_rnaseq.csv"
      all_cell_unmerge_file = "all_cells_unmerge_rnaseq.csv"
      hmec_merge_file = "hmec_rnaseq.csv"
      hmec_unmerge_file = "hmec_unmerge_rnaseq.csv"
      k562_merge_file = "k562_rnaseq.csv"
      k562_unmerge_file = "k562_unmerge_rnaseq.csv"
      mcf7_merge_file = "mcf7_rnaseq.csv"
      mcf7_unmerge_file = "mcf7_unmerge_rnaseq.csv"
      
      # construire la ligne de commande pour le calcul du rnaseq
      # Rscript rnaseq.R chr4 84100000 84600000 PENNY ../merge_rnaseq_cells
      command_line = "Rscript rnaseq.R "
      
      # position
      command_line = paste0(command_line, input$chr, " ", input$position_min,
                            " ", input$position_max)
      
      # unique ID
      command_line = paste0(command_line, " ", tempid)
      
      # si merge file : 1 seul commande
      if(input$merge_rnaseq_file) {
        if (input$merge_rnaseq_experiment){
          command_line = paste0(command_line, "_all ", all_cell_merge_file)
        } else {
          command_line = paste0(command_line, "_all ", all_cell_unmerge_file)
        }
      } else { # sinon 3 commandes, 1 par cell type
        root_command_line <- command_line
        if (input$merge_rnaseq_experiment){
          command_line = paste0(command_line, "_hmec ", hmec_merge_file)
          command_line = paste0(command_line, ";",root_command_line, "_k562 ", k562_merge_file)
          command_line = paste0(command_line, ";",root_command_line, "_mcf7 ", mcf7_merge_file)
        } else {
          command_line = paste0(command_line, "_hmec ", hmec_unmerge_file)
          command_line = paste0(command_line, ";",root_command_line, "_k562 ", k562_unmerge_file)
          command_line = paste0(command_line, ";",root_command_line, "_mcf7 ", mcf7_unmerge_file)
        }
      }
      
      if(isOpen(todornaseqcon,"w")) {
        writeLines(command_line, todornaseqcon)
        output$addrnaseq_msg <- renderText({
          return("Adding RNA-Seq analysis")})
        updateButton(session, "runRNASeq", disabled = FALSE)
      } else {
        output$addrnaseq_msg <- renderText({
          return("Fail adding RNA-Seq analysis")})
      }
      
      
    })
  })
  
  
  observe({
    if(input$runRNASeq > 0) {
      updateButton(session, "addRNASeq", disabled = TRUE)
      updateButton(session, "runRNASeq", disabled = TRUE)
      
      # ENVOI DE MAIL DE CONFIRMATION DE SOUMISSION + MAIL POUR MOI
      # Mail de confirmation de soummision
      
      from <- "no-reply@ShinySNP"
      to <- OPERATOR
      subject <- "RNA-Seq Submission"
      msg <- paste0(USER," has just submitted a RNA-Seq request with id:",tempid,".")
      sendmail(from, to, subject, msg)
      
      from <- "no-reply@ShinySNP"
      to <- USERMAIL
      subject <- "RNA-Seq Submission"
      msg <- paste0("Your request has just been submitted with id:",tempid,".")
      sendmail(from, to, subject, msg)
      
      output$runrnaseq_msg <- renderText({
        return("Sending RNA-Seq analysis request")})
    }
  })
  
  observe({
    
    if(input$runCHIPSeq == 0) {
      return()
    }
    
    isolate ({
      updateButton(session, "runCHIPSeq", disabled = TRUE)
      
      # working files
      hmec_file = "hmec_chipseq.csv"
      k562_file = "k562_chipseq.csv"
      mcf7_file = "mcf7_chipseq.csv"
      
      command_line = paste0("Rscript chipseq.R ", input$chr, " ", input$position_min,
                            " ", input$position_max, " ", tempid)
      
      root_command_line <- command_line
      command_line = paste0(command_line, "_hmec ", hmec_file)
      command_line = paste0(command_line, ";",root_command_line, "_k562 ", k562_file)
      command_line = paste0(command_line, ";",root_command_line, "_mcf7 ", mcf7_file)
      
      # ENVOI DE MAIL DE CONFIRMATION DE SOUMISSION + MAIL POUR MOI
      if(isOpen(todornaseqcon, "w")) {
        writeLines(command_line, todochipseqcon)
        output$runchipseq_msg <- renderText({
          return("Sending ChIP-Seq analysis request")})
        
        # ENVOI DE MAIL DE CONFIRMATION DE SOUMISSION + MAIL POUR MOI
        # Mail de confirmation de soummision
        
        from <- "no-reply@ShinySNP"
        to <- OPERATOR
        subject <- "ChIP-Seq Submission"
        msg <- paste0(USER," has just submitted a ChIP-Seq request with id:",tempid,".")
        sendmail(from, to, subject, msg)
        
        from <- "no-reply@ShinySNP"
        to <- USERMAIL
        subject <- "ChIP-Seq Submission"
        msg <- paste0("Your request has just been submitted with id:",tempid,".")
        sendmail(from, to, subject, msg)
        
      } else {
        output$runchipseq_msg <- renderText({
          return("Fail adding ChIP-Seq analysis")})
      }
    })
  })
  
  
  observe({
    if(input$reset > 0) {
      print("RESET ALL VALUE AND CLEAN ALL FILE CONNECTION")
    }
  })
  
  observe({
    if(input$end > 0) {
      updateButton(session, "end", disabled = TRUE)
      
      if(!is.null(todornaseqcon) && isOpen(todornaseqcon, "w")){close(todornaseqcon)}
      
      if(!is.null(todochipseqcon) && isOpen(todochipseqcon, "w")){close(todochipseqcon)}
      
      if(!is.null(logcon) && isOpen(logcon, "w")){
        writeLines(c(format(Sys.time(), "%a %b %d %X %Y")), logcon)
        close(logcon)
      }
      
      H5close()
      # TROUVER UNE MANIERE MOINS LAIDE DE FERMER L'APPLI
      
      stopApp(42)
    }
  })
  
})

