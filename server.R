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
          updateButton(session, "runRNASeq", disabled = FALSE)
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
        updateSelectInput(session, "chr", selected = "chr12")
        updateNumericInput(session, "position_min", value = 27950000)
        updateNumericInput(session, "position_max", value = 28735000)
        updateNumericInput(session, "snp_position_min", value = 28155080)
        updateNumericInput(session, "snp_position_max", value = 28155080)
        updateNumericInput(session, "hgstart", value = 28111017)
        updateNumericInput(session, "hgend", value = 28127138)
      })}
  })
  
  
  output$highlight <- renderUI({
    list(
      h4("Select regions to highlight"),
      fluidRow(
        column(width = 6,
               numericInput(inputId = "hgstart", 
                            label = "start", 
                            value = NA)),
        column(width = 6,
               numericInput(inputId = "hgend", 
                            label = "end", 
                            value = NA))
      ),
      br(),
      div(style="form-group shiny-input-container",
          selectInput(inputId = "hgcolor", label = "Choose a color", 
                      choices = list("blue" = "steelblue",
                                     "green" = "chartreuse4",
                                     "orange" = "darkorange",
                                     "violet" = "darkviolet",
                                     "pink" = "deeppink",
                                     "red" = "red3",
                                     "yellow" = "gold"))
      ),
      bsButton(inputId = "addhg", label = "Add new Highlight Zone"),
      br(),
      br())
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
        
        msg = paste0("File format error : waiting for a tab delimited file containing ",
                     "at least the following 4 columns : id, start, end, metadata")
        
        stop("File format error : waiting for a tab delimited file containing ",
             "at least the following 4 columns : id, start, end, metadata")
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
      closeAlert(session, alertId = "alert_chr")
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, anchorId = "chr_i", alertId = "alert_chr", title = "Error",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return(waiting_plot("Waiting for your request..."))
        
      } else {
        if(input$position_min >= input$position_max) {
          error_msg <- "Error while setting the region to study : stop position has to be greater than start postion"
          iserror <- TRUE
          createAlert(session, anchorId = "chr_i", alertId = "alert_chr", title = "Error",
                      content = error_msg, append = TRUE, dismiss = TRUE)
          run <<- FALSE
          return(waiting_plot("Waiting for your request..."))
        } else {
          closeAlert(session, alertId = "alert_chr")
          run <<- TRUE
          reset <<- FALSE
        }
      }
    })
  })
  
  drawPlot1 <- reactive({
    
    if (input$run == 0)
      return(waiting_plot("Waiting for your request..."))
    
    if (! "interne" %in% input$my.dataset)
      return(waiting_plot("Waiting for your request..."))
    
    input$reset
    
    if (!run || reset)
      return(waiting_plot("Waiting for your request..."))
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return(waiting_plot("Waiting for your request..."))
        
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
      
      my.data <- load3DData(current_chr = input$chr, my.dataset = "interne")
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
  
  drawPlot1b <- reactive({
    
    if (input$run == 0)
      return(waiting_plot("Waiting for your request..."))
    
    if (! "externe" %in% input$my.dataset)
      return(waiting_plot("Waiting for your request..."))
    
    input$reset
    
    if (!run || reset)
      return(waiting_plot("Waiting for your request..."))
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return(waiting_plot("Waiting for your request..."))
        
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
      
      my.data <- load3DData(current_chr = input$chr, my.dataset = "externe")
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
  
  output$plot1b <- renderPlot({
    drawPlot1b()
  })
  
  
  # affichage de boutons download
  output$download_plot1b <- renderUI ({
    
    if (input$run == 0)
      return()
    
    input$reset
    
    if (!run || reset)
      return()
    
    list(
      br(),
      fluidRow(
        column(width = 6,
               downloadButton(outputId = "download_plot1b_png", label = "Download PNG")),
        column(width = 6,
               downloadButton(outputId = "download_plot1b_pdf", label = "Download PDF"))
      )
    )
  })
  
  
  
  drawPlot23 <- reactive({
    if (input$run == 0)
      return(
        list(plot2 = waiting_plot("Waiting for your request..."), 
             plot3 = waiting_plot("Waiting for your request...")))
    
    input$reset
    
    if (!run || reset)
      return(
        list(plot2 = waiting_plot("Waiting for your request..."), 
             plot3 = waiting_plot("Waiting for your request...")))
    
    isolate ({
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return(waiting_plot("Waiting for your request..."))
        
      } else {
        closeAlert(session, "exampleAlert0")
      }
      
      
      # go ! 
      
      updateButton(session, "reset", disabled = FALSE)
      updateButton(session, "run", disabled = TRUE)
      updateButton(session, "addsnp", disabled = TRUE)
      updateButton(session, "addhg", disabled = TRUE)
      
      print("RUN Plot23")
      
      tryCatch ({
        close(con = snpcon)
        close(hgcon)
      },error = function(e) {
        print(e)
      })
      
      my.data <- load1DData(current_chr = input$chr)
      my.data2 <- load2DData(current_chr = input$chr)
      
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      
      my.ranges <- get1DDataOverview(my.data = my.data, current_range = current_range)
      my.ranges2 <- get2DDataOverview(my.data = my.data2, current_range = current_range)
      
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
      
      impets_tracks <- drawIMPET(ranges_list = my.ranges2, 
                                 current_range = current_range,
                                 highlight_ranges = my.hg.ranges)
      
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      lncrna_figures = drawLNCRNAFigures(my.df = my.data$lncrna, 
                                         current_range = current_range, 
                                         highlight_ranges = my.hg.ranges) 
      
      my.tracks = c(segments_tracks, impets_tracks)
      
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
  
  output$download_plot1b_png <- downloadHandler(
    filename = function() {
      "shinysnp_3D_2.png"
    },
    content = function(file) {
      ggsave(file,drawPlot1b())
    }
  )
  
  output$download_plot1b_pdf <- downloadHandler(
    filename = function() {
      "shinysnp_3D_2.pdf"
    },
    content = function(file) {
      ggsave(file,drawPlot1b())
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
  
  output$download_plot4_png <- downloadHandler(
    filename = function() {
      "shinysnp_rnaseq.png"
    },
    content = function(file) {
      ggsave(file,drawPlot4())
    }
  )
  
  output$download_plot4_pdf <- downloadHandler(
    filename = function() {
      "shinysnp_rnaseq.pdf"
    },
    content = function(file) {
      ggsave(file,drawPlot4())
    }
  )
  
  output$download_plot5_png <- downloadHandler(
    filename = function() {
      "shinysnp_chipseq.png"
    },
    content = function(file) {
      ggsave(file,drawPlot5())
    }
  )
  
  output$download_plot5_pdf <- downloadHandler(
    filename = function() {
      "shinysnp_chipseq.pdf"
    },
    content = function(file) {
      ggsave(file,drawPlot5())
    }
  )
  
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
        msg <- "This process may take several minutes, please be patient and DO NOT refresh the page"
        
        #updateButton(session, "runRNASeq", disabled = TRUE)
        updateButton(session, "reset", disabled = FALSE)
        reset <<- FALSE
        
        tryCatch ({
          close(con = snpcon)
          close(con = hgcon)
        },error = function(e) {
          print(e)
        })
        
      }
      
      closeAlert(session, alertId = "runrnaseq_msg")
      createAlert(session, anchorId = "runrnaseq_msg_i", alertId = "runrnaseq_msg",
                  content = msg, append = FALSE)
    })
  })
  
  drawPlot4 <- reactive({
    
    if (input$runRNASeq == 0)
      return(waiting_plot("Waiting for your request..."))
    
    input$reset
    
    if (!rnaseq || reset)
      return(waiting_plot("Waiting for your request..."))
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return(waiting_plot("Waiting for your request..."))
        
      } else {
        closeAlert(session, "exampleAlert0")
      }
      
      
      # go ! 
      
      updateButton(session, "reset", disabled = FALSE)
      updateButton(session, "run", disabled = TRUE)
      updateButton(session, "addsnp", disabled = TRUE)
      updateButton(session, "addhg", disabled = TRUE)
      
      print("RUN Plot4")
      tryCatch ({
        close(con = snpcon)
        close(hgcon)
      },error = function(e) {
        print(e)
      })
      
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      
      all_cell_merge_file = "all_cells_rnaseq.csv"
      all_cell_unmerge_file = "all_cells_unmerge_rnaseq.csv"
      
      if (input$merge_rnaseq_experiment){
        #localhost
        if(session$clientData$url_hostname == "127.0.0.1") {
          file_list = "/data/hmec_merge_rnaseq.csv"
        } else {
          file_list = paste0("materials/", all_cell_merge_file)
        }        
      } else {
        if(session$clientData$url_hostname == "127.0.0.1") {
          file_list = "/data/hmec_unmerge_rnaseq.csv"
        } else {
          file_list = paste0("materials/", all_cell_unmerge_file)
        }
      }
      
      
      highlight_file = NULL # pas de hg dans ce type de graphe
      t0 = Sys.time()
      rnaseq_tracks <- drawRNASEQ(file_list, highlight_file, current_range)
      print(Sys.time()-t0)
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      ### READ SNP FILE THEN CONVERT IN DATAFRAME
      snpsdf = read.table(snpsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      print(snpsdf)
      
      if(nrow(snpsdf) > 0 ){
        snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
        my.tracks = c(rnaseq_tracks, snp_track, annot_track)
      } else {
        my.tracks = c(rnaseq_tracks, annot_track)
      }
      
      warning("DONE calculate the tracks",call. = FALSE)
      
      tracks(my.tracks) + xlim(current_range)
    })
    
  })
  
  output$plot4 <- renderPlot({
    drawPlot4()
  })
  
  
  # affichage de boutons download
  output$download_plot4 <- renderUI ({
    
    if (input$runRNASeq == 0)
      return()
    
    input$reset
    
    if (!rnaseq || reset)
      return()
    
    list(
      br(),
      fluidRow(
        column(width = 6,
               downloadButton(outputId = "download_plot4_png", label = "Download PNG")),
        column(width = 6,
               downloadButton(outputId = "download_plot4_pdf", label = "Download PDF"))
      )
    )
  })
  
  observe({
    if(input$runCHIPSeq == 0)
      return()
    
    isolate({
      chipseq <<- TRUE
    })
  })
  
  
  observe({
    
    if(input$runCHIPSeq == 0)
      return()
    
    if(!chipseq)
      return()
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        chipseq <<- FALSE
        
        
      } else {
        msg <- "This process may take several minutes, please be patient and DO NOT refresh the page"
        
        #updateButton(session, "runCHIPSeq", disabled = TRUE)
        updateButton(session, "reset", disabled = FALSE)
        reset <<- FALSE
        
        tryCatch ({
          close(con = snpcon)
          close(con = hgcon)
        },error = function(e) {
          print(e)
        })
        
      }
      
      closeAlert(session, alertId = "runchipseq_msg")
      createAlert(session, anchorId = "runchipseq_msg_i", alertId = "runchipseq_msg",
                  content = msg, append = FALSE)
    })
  })
  
  
  drawPlot5 <- reactive({
    
    if (input$runCHIPSeq == 0)
      return(waiting_plot("Waiting for your request..."))
    
    input$reset
    
    if (!chipseq || reset)
      return(waiting_plot("Waiting for your request..."))
    
    isolate ({
      
      if(is.na(input$position_min) || is.na(input$position_max)) {
        error_msg <- "You need to define at least the region to study (chromosome, start and stop position)"
        iserror <- TRUE
        createAlert(session, "alert0", "exampleAlert0", title = "Warning",
                    content = error_msg, append = TRUE, dismiss = TRUE)
        run <<- FALSE
        return(waiting_plot("Waiting for your request..."))
        
      } else {
        closeAlert(session, "exampleAlert0")
      }
      
      
      # go ! 
      
      updateButton(session, "reset", disabled = FALSE)
      updateButton(session, "run", disabled = TRUE)
      updateButton(session, "addsnp", disabled = TRUE)
      updateButton(session, "addhg", disabled = TRUE)
      
      print("RUN Plot5")
      tryCatch ({
        close(con = snpcon)
        close(hgcon)
      },error = function(e) {
        print(e)
      })
      
      current_range <- setStudyRange(current_chr = input$chr, 
                                     current_start = input$position_min, 
                                     current_stop = input$position_max)
      
      #localhost
      if(session$clientData$url_hostname == "127.0.0.1") {
        file_list = paste0("/data/", input$chipseq_cell, "_chipseq.csv")
      } else {
        file_list = paste0("materials/", input$chipseq_cell, "_chipseq.csv")
      }
      
      
      highlight_file = NULL # pas de hg dans ce type de graphe
      t0 = Sys.time()
      rnaseq_tracks <- drawCHIPSEQ(file_list, highlight_file, current_range)
      print(Sys.time()-t0)
      annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
      
      ### READ SNP FILE THEN CONVERT IN DATAFRAME
      snpsdf = read.table(snpsfile, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
      print(snpsdf)
      
      if(nrow(snpsdf) > 0 ){
        snp_track <- drawSNP(current_range = current_range, snps_df = snpsdf, label = "SNPs")
        my.tracks = c(rnaseq_tracks, snp_track, annot_track)
      } else {
        my.tracks = c(rnaseq_tracks, annot_track)
      }
      
      warning("DONE calculate the tracks",call. = FALSE)
      
      tracks(my.tracks) + xlim(current_range)
    })
    
  })
  
  output$plot5 <- renderPlot({
    drawPlot5()
  })
  
  
  # affichage de boutons download
  output$download_plot5 <- renderUI ({
    
    if (input$runCHIPSeq == 0)
      return()
    
    input$reset
    
    if (!chipseq || reset)
      return()
    
    list(
      br(),
      fluidRow(
        column(width = 6,
               downloadButton(outputId = "download_plot5_png", label = "Download PNG")),
        column(width = 6,
               downloadButton(outputId = "download_plot5_pdf", label = "Download PDF"))
      )
    )
  })
  
  
  
  observeEvent(input$reset, {
    # renvoyer sur le panel principal
    session$sendCustomMessage("myCallbackHandler", "1")
    
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
    updateButton(session, "runCHIPSeq", disabled = FALSE)
    updateButton(session, "runRNASeq", disabled = FALSE)
    
    # RESET GLOBALS
    run <<- FALSE
    chipseq <<- FALSE
    rnaseq <<- FALSE
    
    # CLOSING ALL ALERT BOX
    closeAlert(session, alertId = "alert_chr")
    closeAlert(session, alertId = "addhg_msg")
    closeAlert(session, alertId = "addsnp_msg")
    closeAlert(session, alertId = "loadsnp_msg")
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

