# SHINYSNP ui.R
library(shinyBS)

# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)

shinyUI(fluidPage(
  titlePanel("ShinySNP"),
  
  sidebarLayout(
    sidebarPanel(
      helpText(p("1/ Select region"),
               p("2/ Add SNPs"),
               p("3/ Set advanced options (optional)"),
               p("4/ Search for experiments"),
               p("5/ Draw !")),
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
        bsAlert("alert1"),
        bsAlert("alert2"),
        bsAlert("alert3"),
        textOutput("text1"),
        actionButton(inputId = "endAnalysis",label = "End Analysis")

      )
    )
  ))
  