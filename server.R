# SHINYSNP server.R

# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)

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
  
  output$snplist <- renderText({
    snpstart <- input$snp_position_min
    snpend <- input$snp_position_max
    
#     if(is.na(num1) | is.na(num2)) {
#       createAlert(session, "alert", "exampleAlert", title = "Oops",
#                   content = "Both inputs should be numeric.", append = FALSE)
#     } else if(num2 == 0) {
#       createAlert(session, "alert", "exampleAlert", title = "Oops",
#                   content = "You cannot divide by 0.", append = FALSE)
#     } else {
#       closeAlert(session, "exampleAlert")
#       return(num1/num2)
#     }
    return('')
  })
  
})