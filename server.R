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

##### GLOBALS

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

# buttons
run = FALSE
reset = FALSE
rnaseq = FALSE
chipseq = FALSE

# create empty plot
waiting_plot <- function(msg) {
  df = data.frame(x=c(1), 
                  y=c(1), 
                  name = c(msg))
  g = ggplot(data=df, mapping=aes(x=x, y=y)) +
    geom_blank() + ylab("") + xlab("") + 
    geom_text(aes(x = x, y = y, label=name), size = 7, color = "darkgreen") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())
  g
}


setup_files <- function() {
  
  logfile <<- tempfile(pattern = "log_", tmpdir = "logs", fileext = "")
  tempid <<- strsplit(x = logfile, split = "logs/log_", fixed = TRUE)[[1]][2]
  snpsfile <<- paste0("todo/snp_", tempid)
  hgsfile <<- paste0("todo/hg_", tempid)
  todornaseqfile <<- paste0("todo/rnaseq_", tempid)
  todochipseqfile <<- paste0("todo/chipseq_", tempid)
  
  if(!file.exists(logfile)) {
    file.create(logfile)
    print("create log file")
    logcon <<- file(logfile,open="w")
    print("open connection to log file")
    writeLines(c(USERNAME), logcon)
    writeLines(c(format(Sys.time(), "%a %b %d %X %Y")), logcon)
    print("write in log file")
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
  },error = function(e) {
    print(e)
  })
  
  tryCatch ({
    close(hgcon)  
  },error = function(e) {
    print(e)
  })
  
  
  tryCatch ({
    close(todornaseqcon)
  },error = function(e) {
    print(e)
  })
  
  tryCatch ({
    close(todochipseqcon)
  },error = function(e) {
    print(e)
  })
  
  tryCatch ({
    if(isOpen(logcon, "w")){
      writeLines(c(format(Sys.time(), "%a %b %d %X %Y")), logcon)
      print("write in log file before closing")
      close(logcon)
      print("close log file")
    }
    H5close()    
  },error = function(e) {
    print(e)
  })
  
  
  tryCatch ({
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
        h3("Please sign in")
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
          updateButton(session, "reset", disabled = TRUE)
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
    print("addSNPs")
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
    
    if(is.na(selected_chr_min) || is.na(selected_chr_max)) {
      error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
      iserror <- TRUE
      createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                  content = error_msg, append = TRUE, dismiss = TRUE)
    } else {
      closeAlert(session, "exampleAlert0")
    }
    
    
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
    
    if(!is.na(selected_chr_min) && !is.na(selected_chr_max)){
      if(snpstart >= selected_chr_min && snpstart <= selected_chr_max &&
           snpend >= selected_chr_min && snpend <= selected_chr_max) {
        closeAlert(session, "exampleAlert3")          
      } else {
        error_msg <- "The variant must be inclued in the selected region"
        iserror <- TRUE
        createAlert(session, "alert3", "exampleAlert3", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
      }
    }
    
    if(!iserror) {
      writeLines(text = paste(input$snp_label,snpstart,snpend,"n",sep="\t"),
                 con = snpcon)
      updateButton(session, "loadfile", disabled = TRUE)
      
      
      msg = paste0("Adding SNP ",input$snp_label," position : [", snpstart,
                   "-",snpend, "]")
    }
    else
    {
      msg = "Fail while Adding SNP "
    }
    
    closeAlert(session, alertId = "addsnp_msg")
    createAlert(session, anchorId = "addsnp_msg_i", alertId = "addsnp_msg",
                content = msg, append = FALSE)
  })
  
  observe({
    addSNPs()
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
        
        
        msg = paste0("Importing following snps : ", paste(imported_snp$id, 
                                                          collapse = ", "))
        
      }
    },
    error = function(e) {
      print(e)
      msg = paste0("Error while trying to import snps file : ", e)
    })
    
    closeAlert(session, alertId = "loadsnp_msg")
    createAlert(session, anchorId = "loadsnp_msg_i", alertId = "loadsnp_msg",
                content = msg, append = FALSE)
  })
  
  
  observe({
    loadSNPs()
  })
  
  
  addHGs <- eventReactive(input$addhg,{
    print("addHGs")
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
    
    if(is.na(selected_chr_min) || is.na(selected_chr_max)) {
      error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
      iserror <- TRUE
      createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                  content = error_msg, append = TRUE, dismiss = TRUE)
    } else {
      closeAlert(session, "exampleAlert0")
    }
    
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
    
    if(!is.na(selected_chr_min) && !is.na(selected_chr_max)){
      if(hgstart >= selected_chr_min && hgstart <= selected_chr_max &&
           hgend >= selected_chr_min && hgend <= selected_chr_max) {
        closeAlert(session, "exampleAlert6")          
      } else {
        error_msg <- "The HG Zone must be inclued in the selected region"
        iserror <- TRUE
        createAlert(session, "alert6", "exampleAlert6", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
      }
    }
    
    if(!iserror) {
      writeLines(text = paste(hgstart,hgend,input$hgcolor,sep="\t"),
                 con = hgcon)
      msg = paste0("Adding Highlight zone at position [", hgstart,"-",hgend, "] with color : ",input$hgcolor)
      
    } else {
      msg = "Fail while adding Highlight zone"
    }
    closeAlert(session, alertId = "addhg_msg")
    createAlert(session, anchorId = "addhg_msg_i", alertId = "addhg_msg",
                content = msg, append = FALSE)
  })
  
  observe({
    addHGs()
  })
  
  
  observe({
    
    if(input$run == 0)
      return()
    
    isolate({
      run <<- TRUE
      reset <<- FALSE
    })
  })
  
  drawPlot1 <- reactive({
    
    if (input$run == 0)
      return(waiting_plot("Waiting for your request..."))
    
    input$reset
    
    if (!run || reset)
      return()
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return()
        
      } else {
        closeAlert(session, "exampleAlert0")
      }
      
      
      # go ! 
      
      updateButton(session, "reset", disabled = FALSE)
      updateButton(session, "run", disabled = TRUE)
      updateButton(session, "addsnp", disabled = TRUE)
      updateButton(session, "addhg", disabled = TRUE)
      
      print("RUN Plot1")
      tryCatch ({
        close(con = snpcon)
        close(hgcon)
      },error = function(e) {
        print(e)
      })
      
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
  
  
  # affichage de boutons download
  output$download_plot1 <- renderUI ({
    
    if (input$run == 0)
      return()
    
    input$reset
    
    if (!run || reset)
      return()
    
    list(
      br(),
      fluidRow(
        column(width = 6,
               downloadButton(outputId = "download_plot1_png", label = "Download PNG")),
        column(width = 6,
               downloadButton(outputId = "download_plot1_pdf", label = "Download PDF"))
      )
    )
  })
  
  drawPlot23 <- reactive({
    if (input$run == 0)
      return(
        list(plot2 = waiting_plot("Waiting for your request..."), 
             waiting_plot("Waiting for your request...")))
    
    input$reset
    
    if (!run || reset)
      return()
    
    isolate ({
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return()
        
      } else {
        closeAlert(session, "exampleAlert0")
      }
      
      
      # go ! 
      
      updateButton(session, "reset", disabled = FALSE)
      updateButton(session, "run", disabled = TRUE)
      updateButton(session, "addsnp", disabled = TRUE)
      updateButton(session, "addhg", disabled = TRUE)
      
      print("RUN Plot23")
      
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
  
  output$download_plot2 <- renderUI ({
    
    if (input$run == 0)
      return()
    
    input$reset
    
    if (!run || reset)
      return()
    
    list(
      br(),
      fluidRow(
        column(width = 6,
               downloadButton(outputId = "download_plot2_png", label = "Download PNG")),
        column(width = 6,
               downloadButton(outputId = "download_plot2_pdf", label = "Download PDF"))
      )
    )
  })
  
  
  output$plot3 <- renderPlot({
    drawPlot23()$plot3
  })
  
  output$download_plot3 <- renderUI ({
    
    if (input$run == 0)
      return()
    
    input$reset
    
    if (!run || reset)
      return()
    
    list(
      br(),
      fluidRow(
        column(width = 6,
               downloadButton(outputId = "download_plot3_png", label = "Download PNG")),
        column(width = 6,
               downloadButton(outputId = "download_plot3_pdf", label = "Download PDF"))
      )
    )
  })
  
  
  output$download_plot1_png <- downloadHandler(
    filename = function() {
      "shinysnp_3D.png"
    },
    content = function(file) {
      ggsave(file,drawPlot1())
    }
  )
  
  output$download_plot1_pdf <- downloadHandler(
    filename = function() {
      "shinysnp_3D.pdf"
    },
    content = function(file) {
      ggsave(file,drawPlot1())
    }
  )
  
  output$download_plot2_png <- downloadHandler(
    filename = function() {
      "shinysnp_1D.png"
    },
    content = function(file) {
      ggsave(file,drawPlot23()$plot2)
    }
  )
  
  output$download_plot3_png <- downloadHandler(
    filename = function() {
      "shinysnp_hist.png"
    },
    content = function(file) {
      ggsave(file,drawPlot23()$plot3)
    }
  )
  
  output$download_plot2_pdf <- downloadHandler(
    filename = function() {
      "shinysnp_1D.pdf"
    },
    content = function(file) {
      ggsave(file,drawPlot23()$plot2)
    }
  )
  
  output$download_plot3_pdf <- downloadHandler(
    filename = function() {
      "shinysnp_hist.pdf"
    },
    content = function(file) {
      ggsave(file,drawPlot23()$plot3)
    }
  )
  
  observe ({
    
    if(input$addRNASeq == 0) {
      return()
    }
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        msg <- "You need to define at least the region to study (chromosome, start and stop position)"
      } else {
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
        # Rscript src/rnaseq.R chr4 84100000 84600000 PENNY ../merge_rnaseq_cells
        command_line = "Rscript src/rnaseq.R "
        
        # position
        command_line = paste0(command_line, input$chr, " ", input$position_min,
                              " ", input$position_max)
        
        # unique ID
        command_line = paste0(command_line, " ", tempid)
        
        # si merge file : 1 seul commande
        if(input$merge_rnaseq_file) {
          if (input$merge_rnaseq_experiment){
            command_line = paste0(command_line, "_all ", all_cell_merge_file, " RNA-Seq_ALLCELLS todo/snp_", tempid," todo/hg_", tempid)
          } else {
            command_line = paste0(command_line, "_all ", all_cell_unmerge_file, " RNA-Seq_ALLCELLS todo/snp_", tempid," todo/hg_", tempid)
          }
        } else { # sinon 3 commandes, 1 par cell type
          root_command_line <- command_line
          if (input$merge_rnaseq_experiment){
            command_line = paste0(command_line, "_hmec ", hmec_merge_file," RNA-Seq_HMEC todo/snp_", tempid," todo/hg_", tempid)
            command_line = paste0(command_line, ";",root_command_line, "_k562 ", k562_merge_file," RNA-Seq_K562 todo/snp_", tempid," todo/hg_", tempid)
            command_line = paste0(command_line, ";",root_command_line, "_mcf7 ", mcf7_merge_file," RNA-Seq_MCF7 todo/snp_", tempid," todo/hg_", tempid)
          } else {
            command_line = paste0(command_line, "_hmec ", hmec_unmerge_file," RNA-Seq_HMEC todo/snp_", tempid," todo/hg_", tempid)
            command_line = paste0(command_line, ";",root_command_line, "_k562 ", k562_unmerge_file," RNA-Seq_K562 todo/snp_", tempid," todo/hg_", tempid)
            command_line = paste0(command_line, ";",root_command_line, "_mcf7 ", mcf7_unmerge_file," RNA-Seq_MCF7 todo/snp_", tempid," todo/hg_", tempid)
          }
        }
        
        if(isOpen(todornaseqcon,"w")) {
          writeLines(command_line, todornaseqcon)
          msg = "Adding RNA-Seq analysis"
          updateButton(session, "runRNASeq", disabled = FALSE)
        } else {
          msg = "Fail adding RNA-Seq analysis"
        }
      }
      
      closeAlert(session, alertId = "addrnaseq_msg")
      createAlert(session, anchorId = "addrnaseq_msg_i", alertId = "addrnaseq_msg",
                  content = msg, append = FALSE)
      
    })
  })
  
  
  observe({
    
    if(input$runRNASeq == 0)
      return()
    
    isolate({
      rnaseq <<- TRUE
    })
  })
  
  
  observe({
    if(input$runRNASeq == 0)
      return()
    
    if(!rnaseq)
      return()
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        rnaseq <<- FALSE
        
      } else {
        
        
        
        
        updateButton(session, "addRNASeq", disabled = TRUE)
        updateButton(session, "runRNASeq", disabled = TRUE)
        updateButton(session, "reset", disabled = FALSE)
        
        success = FALSE
        
        tryCatch ({
          close(con = snpcon)
          close(hgcon)
        },error = function(e) {
          print(e)
        })
        
        
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
          attachment1 <- mime_part(todornaseqfile, "command_line.txt")
          attachment2 <- mime_part(snpsfile, paste0("snp_", tempid))
          attachment3 <- mime_part(hgsfile, paste0("hg_", tempid))
          msgWithAttachment <- list(msg,attachment1, attachment2, attachment3)
          sendmail(from = from, to = to, subject = subject, msg = msgWithAttachment)
          
          from <- "no-reply@ShinySNP"
          to <- USERMAIL
          subject <- "RNA-Seq Submission"
          msg <- paste0("Your request has just been submitted with id:",tempid,".")
          sendmail(from, to, subject, msg)
        } else {
          rnaseq <<- FALSE
          msg = "Fail submitting RNA-Seq analysis"
        }
      }
      closeAlert(session, alertId = "runrnaseq_msg")
      createAlert(session, anchorId = "runrnaseq_msg_i", alertId = "runrnaseq_msg",
                  content = msg, append = FALSE)
      
    })
  })
  
  observe({
    if(is.null(input$chipseq_analysis)){
      updateButton(session, inputId = "runCHIPSeq", disabled = TRUE)
    } else {
      updateButton(session, inputId = "runCHIPSeq", disabled = FALSE)
    }
    
  })
  
  
  
  observe({
    if(input$runCHIPSeq == 0)
      return()
    
    isolate({
      chipseq <<- TRUE
    })
  })
  
  
  
  observe({
    
    if(input$runCHIPSeq == 0) {
      return()
    }
    
    if(!chipseq)
      return()
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        chipseq <<- FALSE
      } else {
        updateButton(session, "runCHIPSeq", disabled = TRUE)
        updateButton(session, "reset", disabled = FALSE)
        
        tryCatch ({
          close(con = snpcon)
          close(hgcon)
        },error = function(e) {
          print(e)
        })
        
        success = FALSE
        # working files
        hmec_file = "hmec_chipseq.csv"
        k562_file = "k562_chipseq.csv"
        mcf7_file = "mcf7_chipseq.csv"
        
        command_line = ""
        
        root_command_line <-paste0("Rscript src/chipseq.R ", input$chr, " ", input$position_min,
                                   " ", input$position_max, " ", tempid)
        
        if("HMEC" %in% input$chipseq_analysis)
          command_line = paste0(command_line, ";",root_command_line, "_hmec ", hmec_file, " RNA-Seq_HMEC todo/snp_", tempid," todo/hg_", tempid)
        if("K562" %in% input$chipseq_analysis) 
          command_line = paste0(command_line, ";",root_command_line, "_k562 ", k562_file, " RNA-Seq_K562 todo/snp_", tempid," todo/hg_", tempid)
        if("MCF7" %in% input$chipseq_analysis) 
          command_line = paste0(command_line, ";",root_command_line, "_mcf7 ", mcf7_file, " RNA-Seq_MCF7 todo/snp_", tempid," todo/hg_", tempid)
        
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
          attachment1 <- mime_part(todochipseqfile, "command_line.txt")
          attachment2 <- mime_part(snpsfile, paste0("snp_", tempid))
          attachment3 <- mime_part(hgsfile, paste0("hg_", tempid))
          msgWithAttachment <- list(msg,attachment1, attachment2, attachment3)
          sendmail(from = from, to = to, subject = subject, msg = msgWithAttachment)
          
          from <- "no-reply@ShinySNP"
          to <- USERMAIL
          subject <- "ChIP-Seq Submission"
          msg <- paste0("Your request has just been submitted with id:",tempid,".")
          sendmail(from, to, subject, msg)
        } else {
          chipseq <<- FALSE
          msg = "Fail submitting ChIP-Seq analysis"
        }
      }
      
      closeAlert(session, alertId = "runchipseq_msg")
      createAlert(session, anchorId = "runchipseq_msg_i", alertId = "runchipseq_msg",
                  content = msg, append = FALSE)
    })
  })
  
  observeEvent(input$reset, {
    #Desactiver reset pour eviter double click
    updateButton(session, "reset", disabled = TRUE)
    reset <<- TRUE
    
    print("Im IN")
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
    
    
    # RESET RESULTS
    updateButton(session, "addRNASeq", disabled = FALSE)
    updateButton(session, "runCHIPSeq", disabled = FALSE)
    updateButton(session, "runRNASeq", disabled = TRUE)
    
    # RESET GLOBALS
    run <<- FALSE
    chipseq <<- FALSE
    rnaseq <<- FALSE
    
    # CLOSING ALL ALERT BOX
    closeAlert(session, alertId = "addhg_msg")
    closeAlert(session, alertId = "addsnp_msg")
    closeAlert(session, alertId = "loadsnp_msg")
    closeAlert(session, alertId = "addrnaseq_msg")
    closeAlert(session, alertId = "runrnaseq_msg")
    closeAlert(session, alertId = "runchipseq_msg")
    
    print("Resetting the analysis parameters...")
  })
  
  observe({
    if(input$end > 0) {
      updateButton(session, "end", disabled = TRUE)
      
      if (USER$Logged == TRUE) {
        close_files()
      }
      
    }
  })
  
})

