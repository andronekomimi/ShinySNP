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
  
  
  output$addsnp_msg <- renderText({
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
                      content = error_msg, append = TRUE, dismiss = FALSE)
        } else {
          closeAlert(session, "exampleAlert1")
        }
        
        if(snpend == 0 || snpstart == 0){
          error_msg <- "SNP positions are mandatory"
          iserror <- TRUE
          createAlert(session, "alert2", "exampleAlert2", title = "Warning",
                      content = error_msg, append = TRUE, dismiss = FALSE)
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
                      content = error_msg, append = TRUE, dismiss = FALSE)
        }
        
        if(!iserror) {
          writeLines(text = paste(input$snp_label,snpstart,snpend,"n",sep="\t"),
                     con = snpcon)
          return(paste0("Adding SNP ",input$snp_label," position : [", snpstart,
                        "-",snpend, "]"))
        }
      })
    }
  })
  
  output$highlight <- renderUI({
    list(
      div(style="shiny-date-range-input form-group shiny-input-container",
          tags$label("Higlight zone", `for` = "highlight-zone"), 
          div(style="input-daterange input-group",
              tags$input(id = "hgstart", type = "number", value = "start",class="input-small"),
              tags$span("to",style = "input-group-addon"),
              tags$input(id = "hgend", type = "number", value = "end",class="input-small")
          )
      ),
      br(),
      div(style="shiny-date-range-input form-group shiny-input-container",
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
        
        print(paste("start", hgstart))
        print(paste("end", hgend))
        
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
  
  output$plot <- renderPlot({
    
      if (input$draw == 0)
        return()
      
      isolate ({
        print(input$chr)
        my.data <- loadChrData(current_chr = input$chr)
        current_range <- setStudyRange(current_chr = input$chr, 
                                       current_start = input$position_min, 
                                       current_stop = input$position_max)
        my.ranges <- getDataOverview(my.data = my.data, current_range = current_range)
        print(length(my.ranges$arch))
        
        archs_tracks <- drawArchs(ranges_list = my.ranges$arch, 
                                  current_range = current_range,
                                  highlight_range_list = NULL)
        
        annot_track <- drawAnnotations("Genes",current_range = current_range + 10000)
        #snp_track <- drawSNP(current_range, current_chr, as.numeric(selected_row["Position"]), selected_row["RsID"])
        
        t = c(archs_tracks, annot_track)
        
        print("DONE")
        
        tracks(t) + xlim(current_range)
      })
    
  })

  observe({
    if(input$endAnalysis > 0) {
      updateButton(session, "endAnalysis", disabled = TRUE)
      if(isOpen(snpcon)){close(snpcon)}
      if(isOpen(hgcon)){close(hgcon)}
      if(isOpen(logcon)){
        writeLines(c(format(Sys.time(), "%a %b %d %X %Y")), logcon)
        close(logcon)
      }
      #if(file.exists(snpsfile)) {file.remove(snpsfile)}
    }
  })
  
  
  # subset(sgp_k562_lane13, (sgp_k562_lane13$V2 >= current_start & sgp_k562_lane13$V2 <= current_stop) | (sgp_k562_lane13$V3 >= current_start & sgp_k562_lane13$V3 <= current_stop))
})

