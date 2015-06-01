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
                  selected = "chr1"),
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
      h6("27950000-28735000"),
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
      actionButton(inputId = "addsnp",label = "Add new SNP"),
      hr(),
      h4(" OR "),
      shinyFilesButton('loadsnp', 'Load SNP file', 
                       'Please select a file', FALSE),
      textInput(inputId = "path", label = ""),
      br(),
      br(),
      span(textOutput("addsnp_msg"), style = "color:green"),
      br(),
      br(),
      bsButton(inputId = "draw", label = "Draw", style = "btn btn-primary")
    ),
    
    mainPanel(
      bsAlert("alert1"),
      bsAlert("alert2"),
      bsAlert("alert3"),
      bsAlert("alert4"),
      bsAlert("alert5"),
      bsAlert("alert6"),
      navbarPage(span("Analysis", style = "color:green"),
                 tabPanel("Conformation",
                          #       conditionalPanel(condition = "input.draw > 0 && !input.plot",
                          #                        list(
                          #                          img(src = "page_loader.gif",
                          #                              filetype = "image/gif",
                          #                              alt = "Please wait... I'm processing your query")
                          #                        )
                          #       ),
                          plotOutput("plot1"),
                          br(),
                          br()
                 ),
                 tabPanel("Lncrna & Enhancers",
                          plotOutput("plot2"),
                          br(),
                          br()
                 ),
                 tabPanel("RNA-Seq",
                          h3("RNA-Seq"),
                          helpText("We are working with the following dataset : ",
                                   tags$ul(list(
                                     tags$li("MCF7 : 7 experiments with 2 replicates each 
                                   and 1 experiment with 3 replicates"),
                                     tags$li("K562 : 3 experiments with 2 replicates for each"),
                                     tags$li("HMEC : 1 experiment with 2 replicates")
                                   ))),
                          br(),
                          fluidRow(
                            column(width = 4,
                                   radioButtons(inputId = 'merge_rnaseq_replicate', 
                                                label = "Analysis configuration",
                                                choices = list("One file per cell type" = FALSE, 
                                                               "Merge in a single file" = TRUE),
                                                selected = TRUE)),
                            column(width = 4,
                                   radioButtons(inputId = 'merge_rnaseq_experiment', 
                                                label = "",
                                                choices = list("One track per experiment" = FALSE, 
                                                               "Merge experiments" = TRUE),
                                                selected = TRUE)),
                            column(width = 4,
                                   radioButtons(inputId = 'merge_rnaseq_file', 
                                                label = "",
                                                choices = list("One track per replicate" = FALSE, 
                                                               "Merge replicates" = TRUE),
                                                selected = TRUE))
                          ),
                          br(),
                          br()
                 ),
                 tabPanel("ChIP-Seq",
                          h3("ChIP-Seq"),
                          helpText("We are working with the following dataset : ",
                                   tags$ul(list(
                                     tags$li("MCF7 : 7 experiments with 2 replicates each 
                                   and 1 experiment with 3 replicates"),
                                     tags$li("K562 : 3 experiments with 2 replicates for each"),
                                     tags$li("HMEC : 1 experiment with 2 replicates")
                                   ))),
                          br(),
                          br()
                 )
                 
      ),
      bsButton(inputId = "endAnalysis",label = "End Analysis", style = "btn btn-primary")
      
    )
  )
))
