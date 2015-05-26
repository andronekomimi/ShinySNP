# SHINYSNP ui.R
library(shinyBS)

# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)

shinyUI(fluidPage(
  titlePanel("ShinySNP"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("1/ Select region
                2/ Add SNPs
                3/ Set advanced options (optional)
                4/ Search for experiments
                5/ Draw !"),
      h3("Select region"),
      selectInput("chr", 
                  label = "Choose a chromosome",
                  choices = as.character(chroms$chr),
                  selected = NA),
      fluidRow(
        column(width = 6,
               numericInput(inputId = 'position_min', 
                            label = "start",
                            value = NA)),
        column(width = 6,
               numericInput(inputId = 'position_max', 
                            label = "end",
                            value = NA))),
      br(),
      h3("Add SNPs"),
      textInput(inputId = "snp_label",label = "SNP Label", value = "My SNP"),
      fluidRow(
        column(width = 6,
               numericInput(inputId = 'snp_position_min', 
                            label = "start",
                            value = 0)),
        column(width = 6,
               numericInput(inputId = 'snp_position_max', 
                            label = "end",
                            value = 0))),
      actionButton(inputId = "addsnp",label = "Add new SNP"),
      textOutput("snplist")
    ),
    
    mainPanel(
      bsAlert("alert"),
      textOutput("text1")
    )
  )
))
