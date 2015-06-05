# SHINYSNP ui.R
library(shinyBS)
library(shinyFiles)

# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)

shinyUI(fluidPage(
  titlePanel("Shiny SNP"),
  
  sidebarLayout(
    sidebarPanel(
      #       helpText(p("1/ Select region"),
      #                p("2/ Add SNPs"),
      #                p("3/ Set advanced options (optional)"),
      #                p("4/ Search for experiments"),
      #                p("5/ Draw !")),
      h3("Select region"),
      selectInput("chr", 
                  label = "Choose a chromosome",
                  choices = as.character(chroms$chr),
                  selected = "chr12"),
      fluidRow(
        column(width = 6,
               numericInput(inputId = 'position_min', 
                            label = "start",
                            value = NA)),
        column(width = 6,
               numericInput(inputId = 'position_max', 
                            label = "end",
                            value = NA))
      ),
      bsCollapse(id = "param", open = NULL, multiple = FALSE,
                 bsCollapsePanel('Advanced parameters',uiOutput("highlight"), 
                                 style = "default")
      ),
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
      bsButton(inputId = "addsnp",label = "Add new SNP"),
      hr(),
      h4(" OR "),
      shinyFilesButton('loadsnp', 'Load SNP file', 
                       'Please select a file', FALSE),
      textInput(inputId = "path", label = ""),
      br(),
      br(),
      span(textOutput("addsnp_msg"), style = "color:green"),
      span(textOutput("addsnp_err"), style = "color:red"),
      br(),
      h3("Choose the dataset"),
      radioButtons(inputId = 'my.dataset', 
                   inline = TRUE,
                   label = NULL,
                   choices = list("local" = "interne", 
                                  "4DGenome" = "externe"),
                   selected = "interne"),
      br(),
      fluidRow(
        column(width = 4,
               bsButton(inputId = "run", label = "Run Analysis", style = "btn btn-primary")),
        column(width = 4,
               bsButton(inputId = "reset", label = "New Analysis", style = "btn btn-primary")),
        column(width = 4,
               bsButton(inputId = "end",label = "End Analysis", style = "btn btn-primary"))),
      br(),
      actionButton(inputId = "runEx", label = "Run Example"),
      br()
    ),
    mainPanel(
      bsAlert("alert1"),
      bsAlert("alert2"),
      bsAlert("alert3"),
      bsAlert("alert4"),
      bsAlert("alert5"),
      bsAlert("alert6"),
      navbarPage(span("Analysis", style = "color:green"),
                 tabPanel("3D Conformation",
#                                 conditionalPanel(condition = "input.draw > 0 && output.plot1",
#                                                  list(
#                                                    img(src = "page_loader.gif",
#                                                        filetype = "image/gif",
#                                                        alt = "Please wait... I'm processing your query")
#                                                  )
#                                 ),
                          plotOutput("plot1"),
                          br(),
                          br()
                 ),
                 tabPanel("Lncrna & Enhancers",
                          plotOutput("plot2"),
                          plotOutput("plot3"),
                          br(),
                          br()
                 ),
                 tabPanel("RNA-Seq",
                          h3("RNA-Seq"),
                          helpText("We are working with the following dataset : ",
                                   tags$ul(list(
                                     tags$li("MCF7 : 7 experiments (2 replicates/exp) 
                                   and 1 experiment (3 replicates/exp)"),
                                     tags$li("K562 : 3 experiments (2 replicates/exp)"),
                                     tags$li("HMEC : 1 experiment (2 replicates/exp)")
                                   ))),
                          br(),
                          fluidRow(
                            column(width = 6,
                                   radioButtons(inputId = 'merge_rnaseq_experiment', 
                                                label = "",
                                                choices = list("One track per experiment" = FALSE, 
                                                               "Merge experiments" = TRUE),
                                                selected = TRUE)),
                            column(width = 6,
                                   radioButtons(inputId = 'merge_rnaseq_file', 
                                                label = "",
                                                choices = list("One file per cell type" = FALSE, 
                                                               "One single file" = TRUE),
                                                selected = TRUE))
                          ),
                          br(),
                          fluidRow(
                            column(width = 6,
                                   bsButton(inputId = "addRNASeq", label = "Add RNA-Seq analysis"),
                                   br(), span(textOutput("addrnaseq_msg"), style = "color:green"), br()),
                            column(width = 6,
                                   bsButton(inputId = "runRNASeq",label = "Send analysis request", style = "btn btn-primary", disabled = TRUE),
                                   br(), span(textOutput("runrnaseq_msg"), style = "color:green")))
                 ),
                 tabPanel("ChIP-Seq",
                          h3("ChIP-Seq"),
                          helpText("We are working with the following dataset : ",
                                   tags$ul(list(
                                     tags$li("Targets for MCF7 : CTCF, H3K27ac, H3K36me3, H3K4me3, H3K27me3"),
                                     tags$li("Targets for K562 : H3K36me3, H3K4me2, CTCF, EZH2"),
                                     tags$li("Targets for HMEC : H3K4me3, H2AFZ, H3K9me3, H3K36me3, H3K4me3, 
                                             H3K27ac, H3K4me1, H3K79me2, CTCF, H3K9ac, H3K4me2, 
                                             H4K20me1, H3K27me3, EZH2")
                                   ))),
                          br(),
                          br(),
                          bsButton(inputId = "runCHIPSeq",label = "Send analysis request", style = "btn btn-primary"),
                          br(), br(), span(textOutput("runchipseq_msg"), style = "color:green")
                 )
                 
      )     
    )
  )
))
