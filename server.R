# SHINYSNP server.R

source("data_functions.R")
library(sendmailR)

# EXECUTER 1 FOIS AU LANCEMENT DE LAPPLI
# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)

Logged = FALSE;

## Chargement de la DB d'authenfication
membership_db = "/etc/shiny-apps/ShinySNP_membership.sqlite"
con = RSQLite::dbConnect(RSQLite::SQLite(), membership_db)
members = dbGetQuery(con,'select * from membership')
RSQLite::dbDisconnect(con)

OPERATOR = "audrey.lemacon.1@ulaval.ca"
USERMAIL = ""
USERNAME = ""
USERROLE = ""

# files connections
tempid = 0

logfile = ""
logcon = ""

snpcon = ""
snpsfile = ""

hgcon = ""
hgsfile = ""

todornaseqcon = ""
todornaseqfile = ""

todochipseqcon = ""
todochipseqfile = ""

setup_files <- function() {
  
  logfile <<- tempfile(pattern = "log_", tmpdir = "logs", fileext = "")
  tempid <<- strsplit(x = logfile, split = "logs/log_", fixed = TRUE)[[1]][2]
  snpsfile <<- paste0("logs/snp_", tempid)
  hgsfile <<- paste0("logs/hg_", tempid)
  todornaseqfile <<- paste0("todo/rnaseq_", tempid)
  todochipseqfile <<- paste0("todo/chipseq_", tempid)
  
  if(!file.exists(logfile)) {
    file.create(logfile)
    logcon <<- file(logfile,open="w")
    writeLines(c(USERNAME), logcon)
    writeLines(c(format(Sys.time(), "%a %b %d %X %Y")), logcon)
  }
  
  if(!file.exists(snpsfile)) {
    file.create(snpsfile)
    snpcon <<- file(snpsfile,open="a+b")
    writeLines(paste("id","start","end","metadata",sep = "\t"), snpcon)
  }
  
  if(!file.exists(hgsfile)) {
    file.create(hgsfile)
    hgcon <<- file(hgsfile,open="a+b")
    writeLines(paste("start","end","color",sep = "\t"), hgcon)
  }
  
  if(!file.exists(todornaseqfile)) {
    file.create(todornaseqfile)
    todornaseqcon <<- file(todornaseqfile,open="w")
    writeLines(paste0("#", USERNAME), todornaseqcon)
    writeLines(paste0("#", tempid), todornaseqcon)
  }
  
  if(!file.exists(todochipseqfile)) {
    file.create(todochipseqfile)
    todochipseqcon <<- file(todochipseqfile,open="w")
    writeLines(paste0("#", USERNAME), todochipseqcon)
    writeLines(paste0("#", tempid), todochipseqcon)
  }
  
  print(paste0("Files setup for id:",tempid))
}


close_files <- function() {
  # close all resilient files
  print("closing files")
  tryCatch ({
    close(snpcon)
    close(hgcon)
    close(todornaseqcon)
    close(todochipseqcon)
    if(!is.null(logcon) && isOpen(logcon, "w")){
      writeLines(c(format(Sys.time(), "%a %b %d %X %Y")), logcon)
      close(logcon)
    }
    H5close()    
  },error = function(e) {
    print(e)
  })
}

options(scipen=999) # virer l'annotation scientifique

shinyServer(function(input, output, session) {
  USER <- reactiveValues(Logged = FALSE)
  
  output$uiLogin <- renderUI({
    if (USER$Logged == FALSE) {
      list(
        h3("Please sign in"),
        wellPanel(style = "background-color: #ffffff;",
                  textInput("Username", "User Name:"),
                  passwordInput("Password", "Password:"),
                  br(),
                  actionButton("Login", "Log in"),
                  br(),br(),
                  span(textOutput("pass"), style = "color:red")
        )
      )
    }
  })
  
  outputOptions(output, "uiLogin", suspendWhenHidden=FALSE)
  
  output$pass <- renderText({  
    if (USER$Logged == FALSE) {
      
      if (!is.null(input$Login) && input$Login > 0) {
        authentif_success = TRUE
        
        Username <- isolate(input$Username)
        Password <- isolate(input$Password)
        
        if(!(Username %in% members$user)) {
          authentif_success = FALSE
        } else {
          passw = subset(members, user==Username)$password
          if(passw != Password) {
            authentif_success = FALSE
          }
        }
        
        if (authentif_success) {
          USER$Logged <- TRUE
          USERNAME <<- as.character(subset(members, user==Username)$username)
          USERMAIL <<- as.character(subset(members, user==Username)$email)
          USERROLE <<- as.character(subset(members, user==Username)$role)
          
          setup_files()
          
          # Reactiver ts les boutons
          updateButton(session, "addsnp", disabled = FALSE)
          updateButton(session, "addhg", disabled = FALSE)
          updateButton(session, "run", disabled = FALSE)
          updateButton(session, "reset", disabled = FALSE)
          updateButton(session, "addRNASeq", disabled = FALSE)
          updateButton(session, "runCHIPSeq", disabled = FALSE)
          
          output$runEx <- renderUI({
            if(USERROLE == "admin") {
              return(list(actionButton(inputId = "runEx", label = "Run Example")))
            }
          })
          
          "Authentication succeed!"
          
        } else  {
          "Authentication failed!"
        }
      } 
    }
  })
  
  
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
    if(!is.null(input$runEx)){
      if(input$runEx == 0)
        return()
      
      isolate({
        updateNumericInput(session, "position_min", value = 27950000)
        updateNumericInput(session, "position_max", value = 28735000)
        updateNumericInput(session, "snp_position_min", value = 28155080)
        updateNumericInput(session, "snp_position_max", value = 28155080)
        updateNumericInput(session, "hgstart", value = 28111017)
        updateNumericInput(session, "hgend", value = 28127138)
      })}
  })
  
  
  addSNPs <- eventReactive(input$addsnp, {
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
      updateButton(session, "loadfile", disabled = TRUE)
      
      
      return(paste0("Adding SNP ",input$snp_label," position : [", snpstart,
                    "-",snpend, "]"))
    }
    else
    {
      return("Fail while Adding SNP ")
    }
  })
  
  observe({
    output$addsnp_msg <- renderText({
      addSNPs()
    })
  })
  
  
  loadSNPs <- eventReactive(input$loadfile, {
    inFile <- input$loadfile
    
    if (is.null(inFile))
      return(NULL)
    
    tryCatch ({
      # get data from imported snp file
      imported_snp = read.table(inFile$datapath, header=TRUE, 
                                stringsAsFactors=FALSE, quote = "\"", sep="\t")
      if(is.null(snpcon)) {
        print('ok')
      }
      
      if((!"id" %in% names(imported_snp)) || (!"start" %in% names(imported_snp)) ||
           (!"end" %in% names(imported_snp)) || (!"metadata" %in% names(imported_snp))) {
        stop(paste0("File format error : waiting for a tab delimited file containing ",
                    "at least the following 4 columns : id, start, end, metadata"))
      } else {
        updateButton(session,inputId = "addsnp", disabled = TRUE)
        
        for (i in (1:nrow(imported_snp))) {
          new_snp = imported_snp[i,]
          writeLines(text = paste(new_snp$id,new_snp$start,new_snp$end,
                                  new_snp$metadata,sep="\t"), con = snpcon)
        }
        
        
        return(paste0("Importing following snps : ", paste(imported_snp$id, 
                                                           collapse = ", ")))
        
      }
    },
    error = function(e) {
      print(e)
      return(paste0("Error while trying to import snps file : ", e))
    })
  })
  
  
  observe({
    output$loadsnp_msg <- renderText({
      loadSNPs()
    })
  })
  
  
  addHGs <- eventReactive(input$addhg,{
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
      writeLines(text = paste(hgstart,hgend,input$hgcolor,sep="\t"),
                 con = hgcon)
      return(paste0("Adding Highlight zone at position [", hgstart,
                    "-",hgend, "] with color : ",input$hgcolor))
    }
    
  })
  
  
  observe({
    output$addhg_msg <- renderText({
      if(!is.na(input$position_min)) {     
        print("test !!")
        addHGs()
      }
    })
  })
  
  observe({
    
    if(input$run == 0)
      return()
    
    isolate({
      updateButton(session, "run", disabled = TRUE)
      updateButton(session, "addsnp", disabled = TRUE)
      updateButton(session, "addhg", disabled = TRUE)
      
    })
  })
  
  drawPlot1 <- reactive({
    if (input$run == 0)
      return()
    
    isolate ({
      if(!is.null(snpcon) && isOpen(snpcon, "a+b")){
        close(con = snpcon)
      }
      
      if(!is.null(hgcon) && isOpen(hgcon, "a+b")){
        close(hgcon)
      }
      
      my.data <- load3DData(current_chr = input$chr, my.dataset = input$my.dataset)
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      my.ranges <- get3DDataOverview(my.data = my.data, current_range = current_range)
      
      highlights = read.table(hgsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      print(highlights)
      
      my.hg.ranges = NULL
      
      if(nrow(highlights)){
        my.hg.ranges = GRanges(seqnames = input$chr, 
                               ranges = IRanges(highlights$start,highlights$end),  
                               alpha = rep(0.5,times = nrow(highlights)),
                               color = highlights$color)
      }
      
      archs_tracks <- drawArchs(ranges_list = my.ranges,
                                current_range = current_range,
                                highlight_ranges = my.hg.ranges)
      
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      ### READ SNP FILE THEN CONVERT IN DATAFRAME
      print(snpsfile)
      snpsdf = read.table(snpsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      print(snpsdf)
      
      if(nrow(snpsdf) > 0 ){
        snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
        my.tracks = c(archs_tracks, snp_track, annot_track)
      } else {
        my.tracks = c(archs_tracks, annot_track)
      }
      
      warning("DONE calculate the tracks",call. = FALSE)
      
      tracks(my.tracks) + xlim(current_range)
    })
    
  })
  
  output$plot1 <- renderPlot({
    drawPlot1()
  })
  
  drawPlot23 <- reactive({
    if (input$run == 0)
      return()
    
    isolate ({
      
      my.data <- load1DData(current_chr = input$chr)
      
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      my.ranges <- get1DDataOverview(my.data = my.data, current_range = current_range)
      
      highlights = read.table(hgsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      print(highlights)
      
      my.hg.ranges = NULL
      
      if(nrow(highlights)){
        my.hg.ranges = GRanges(seqnames = input$chr, 
                               ranges = IRanges(highlights$start,highlights$end),  
                               alpha = rep(0.5,times = nrow(highlights)),
                               color = highlights$color)
      }
      
      segments_tracks <- drawSegment(ranges_list = my.ranges, 
                                     current_range = current_range,
                                     highlight_ranges = my.hg.ranges)
      
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      lncrna_figures = drawLNCRNAFigures(my.df = my.data$lncrna, 
                                         current_range = current_range, 
                                         highlight_ranges = my.hg.ranges) 
      
      my.tracks = c(segments_tracks)
      
      if(!is.null(lncrna_figures$lncrna_track)){
        my.tracks = c(my.tracks, lncrna_figures$lncrna_track)
      }
      
      ### READ SNP FILE THEN CONVERT IN DATAFRAME
      snpsdf = read.table(snpsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      
      if(nrow(snpsdf) > 0 ){
        snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
        my.tracks = c(my.tracks, snp_track)
      }
      
      my.tracks = c(my.tracks, annot_track)
      
      warning("DONE calculate the tracks",call. = FALSE)
      
      list(plot2 = tracks(my.tracks) + xlim(current_range), plot3 = tracks(lncrna_figures$lncrna_hist))
    })
  })
  
  output$plot2 <- renderPlot({
    drawPlot23()$plot2
  })
  
  output$plot3 <- renderPlot({
    drawPlot23()$plot3
  })
  
  output$download_plot1 <- downloadHandler(
    filename = function() {
      "shinysnp_3D.png"
    },
    content = function(file) {
      ggsave(file,drawPlot1())
    }
  )
  
  output$download_plot2 <- downloadHandler(
    filename = function() {
      "shinysnp_1D.png"
    },
    content = function(file) {
      ggsave(file,drawPlot23()$plot2)
    }
  )
  
  output$download_plot3 <- downloadHandler(
    filename = function() {
      "shinysnp_hist.png"
    },
    content = function(file) {
      ggsave(file,drawPlot23()$plot3)
    }
  )
  
  output$addsnp_msg <- renderText({
    
    
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
      success = FALSE
      
      if(isOpen(todornaseqcon, "w")) {
        close(todornaseqcon)
        success = TRUE
      } 
      
      
      # ENVOI DE MAIL DE CONFIRMATION DE SOUMISSION + MAIL POUR MOI
      # Mail de confirmation de soummision
      
      if(success){
        output$runrnaseq_msg <- renderText({
          return("Sending RNA-Seq analysis request")})
        
        from <- "no-reply@ShinySNP"
        to <- OPERATOR
        subject <- "RNA-Seq Submission Report"
        msg <- paste0(USERNAME," has just submitted a RNA-Seq request with id:",tempid,".")
        attachment <- mime_part(todornaseqfile, "command_line.txt")
        msgWithAttachment <- list(msg,attachment)
        sendmail(from = from, to = to, subject = subject, msg = msgWithAttachment)
        
        
        from <- "no-reply@ShinySNP"
        to <- USERMAIL
        subject <- "RNA-Seq Submission"
        msg <- paste0("Your request has just been submitted with id:",tempid,".")
        sendmail(from, to, subject, msg)
      } else {
        output$runrnaseq_msg <- renderText({
          return("Fail submitting RNA-Seq analysis")})
      }
      
    }
  })
  
  observe({
    
    if(input$runCHIPSeq == 0) {
      return()
    }
    
    isolate ({
      updateButton(session, "runCHIPSeq", disabled = TRUE)
      success = FALSE
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
      if(isOpen(todochipseqcon, "w")) {
        writeLines(command_line, todochipseqcon)
        close(todochipseqcon)
        success = TRUE
      } 
      
      if(success){
        output$runchipseq_msg <- renderText({
          return("Sending ChIP-Seq analysis request")})
        
        # ENVOI DE MAIL DE CONFIRMATION DE SOUMISSION + MAIL POUR MOI
        # Mail de confirmation de soummision
        
        from <- "no-reply@ShinySNP"
        to <- OPERATOR
        subject <- "ChIP-Seq Submission Report"
        msg <- paste0(USERNAME," has just submitted a ChIP-Seq request with id:",tempid,".")
        attachment <- mime_part(todochipseqfile, "command_line.txt")
        msgWithAttachment <- list(msg,attachment)
        sendmail(from = from, to = to, subject = subject, msg = msgWithAttachment)
        
        from <- "no-reply@ShinySNP"
        to <- USERMAIL
        subject <- "ChIP-Seq Submission"
        msg <- paste0("Your request has just been submitted with id:",tempid,".")
        sendmail(from, to, subject, msg)
      } else {
        output$runchipseq_msg <- renderText({
          return("Fail submitting ChIP-Seq analysis")})
      }
      
      
    })
  })
  
  observeEvent(input$reset, {
    
    # close previous log file
    close_files()
    
    # open new set of files
    setup_files()
    
    # RESET PARAMETERS
    updateNumericInput(session, "position_min", value = NA)
    updateNumericInput(session, "position_max", value = NA)
    updateNumericInput(session, "snp_position_min", value = NA)
    updateNumericInput(session, "snp_position_max", value = NA)
    updateNumericInput(session, "hgstart", value = NA)
    updateNumericInput(session, "hgend", value = NA)
    updateButton(session, "addsnp", disabled = FALSE)
    updateButton(session, "addhg", disabled = FALSE)
    updateButton(session, "run", disabled = FALSE)
    updateButton(session, "loadfile", disabled = FALSE)
    updateButton(session, "reset", disabled = FALSE)
    output$addsnp_msg <- renderText({return("")})
    output$loadsnp_msg <- renderText({return("")})
    output$addhg_msg <- renderText({return("")})
    
    # RESET RESULTS
    
    updateButton(session, "addRNASeq", disabled = FALSE)
    updateButton(session, "runCHIPSeq", disabled = FALSE)
    updateButton(session, "runRNASeq", disabled = TRUE)
    output$addrnaseq_msg <- renderText({return("")})
    output$runrnaseq_msg <- renderText({return("")})
    output$runchipseq_msg <- renderText({return("")})
    
    output$plot1 <- renderPlot({})
    output$plot2 <- renderPlot({})
    output$plot3 <- renderPlot({})
  })
  
  
  observe({
    if(input$end > 0) {
      updateButton(session, "end", disabled = TRUE)
      
      if (USER$Logged == TRUE) {
        close_files()
      }
      # TROUVER UNE MANIERE MOINS LAIDE DE FERMER L'APPLI
      
      stopApp(42)
    }
  })
  
})

