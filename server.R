# SHINYSNP server.R

# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)
con = ""

shinyServer(function(input, output, session) {
  
  reactive ({
    snpsfile <- paste0("logs/", session$clientData$singletons,".snps")
    
    print(snpsfile)
    
#     if(!file.exists(snpsfile)) {
#       
#       file.create(snpsfile)
#       con=file(snpsfile,open="a+b")
#       writeLines("snp_id  start end metadata", con)
#     }
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
  
  
  output$snplist <- renderText({
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
#           writeLines(text = paste(input$snp_label,snpstart,snpend,"n",sep="\t"),
#                      con = con)
          return("LOOOOL")
        }
      })
    }
  })
  
  observe({
    if(input$endAnalysis > 0) {
      #snpsfile <- paste0("logs/", session$clientData$singletons,".snps")
      #if(isOpen(con)){close(con)}
      #if(file.exists(snpsfile)) {file.remove(snpsfile)}
    }
  })
  
  
  
})

