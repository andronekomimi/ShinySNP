# SHINYSNP ui.R
library(shinyBS)

# chrom_list dans les parametres
chroms = read.csv("data/chromosomes.txt", header = TRUE)


footer<-function(){
  tags$div(
    class = "footer",
    tags$div(class = "foot-inner", 
             list(
               "Shiny SNP is an ", 
               tags$a(href="http://bioinformatique.genome.ulaval.ca/recherche/","Arnaud Droit Lab"),
               " project (2015).",
               br(),
               "Functional design and Development ",
               tags$a(href="https://ca.linkedin.com/in/audreylemacon", "Audrey Lemacon ")
             )
             
    )
  )
}


shinyUI(fluidPage(
  tags$head(
    tags$link(rel = "icon", type = "image/x-icon", href = "logo.jpg")
  ),
  titlePanel("Shiny SNP"),
  uiOutput(outputId ="uiLogin"),
  conditionalPanel(condition = "output.uiLogin",
                   list(
                     wellPanel(style = "background-color: #ffffff;",
                               textInput("Username", "User Name:"),
                               passwordInput("Password", "Password:"),
                               br(),
                               actionButton("Login", "Log in"),
                               br(),br(),
                               span(textOutput("pass"), style = "color:red")
                     )
                   )
  ),
  sidebarLayout(
    conditionalPanel(condition = "!output.uiLogin",
                     list(
                       sidebarPanel(

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
                         hr(),
                         h3("Select regions to highlight"),
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
                         br(),
                         bsAlert("addhg_msg_i"),
                         hr(),
                         h3("Add variants"),
                         helpText("Add your variants manually or load a file conaining the variants list"),
                         helpText("File format : ",
                                  tags$ul(
                                    tags$li("Tab-delimited, with 4 named column 'id,start,stop,metadata'"),
                                    tags$li("One variant per line"),
                                    tags$li("The metadata column enables to make variants groups. 
               The variants of the same group will be represent with the same color"))),
                         textInput(inputId = "snp_label",label = "Variant Label", value = "My SNP"),
                         fluidRow(
                           column(width = 6,
                                  numericInput(inputId = 'snp_position_min', 
                                               label = "start",
                                               value = NA)),
                           column(width = 6,
                                  numericInput(inputId = 'snp_position_max', 
                                               label = "end",
                                               value = NA))),
                         bsButton(inputId = "addsnp",label = "Add new SNP"),
                         br(),br(),
                         bsAlert("addsnp_msg_i"),
                         hr(),
                         fileInput('loadfile', 'Choose file to upload',
                                   accept = c(
                                     'text/csv',
                                     'text/comma-separated-values',
                                     'text/tab-separated-values',
                                     'text/plain',
                                     '.csv',
                                     '.tsv',
                                     '.txt',
                                     ''
                                   )
                         ),
                         bsAlert("loadsnp_msg_i"),
                         hr(),
                         h3("Choose the dataset"),
                         radioButtons(inputId = 'my.dataset', 
                                      inline = TRUE,
                                      label = NULL,
                                      choices = list("local" = "interne", 
                                                     "4DGenome" = "externe"),
                                      selected = "interne"),
                         br(),
                         uiOutput(outputId = "runEx"),
                         fluidRow(
                           column(width = 4,
                                  bsButton(inputId = "run", label = "Run Analysis", style = "btn btn-primary", disabled = TRUE)),
                           column(width = 4,
                                  bsButton(inputId = "reset", label = "New Analysis", style = "btn btn-primary", disabled = TRUE)),
                           column(width = 4,
                                  tags$button(type = "button", class="btn btn-default action-button shiny-bound-input",
                                              id = "end", "End Analysis", onClick="history.go(0)")
                           ))
                       ))
    ),
    mainPanel(
      conditionalPanel(condition = "!output.uiLogin",
                       list(
                         navbarPage(span("Analysis", style = "color:green"),
                                    tabPanel("3D Conformation",
                                             bsAlert("alert0"),
                                             bsAlert("alert1"),
                                             bsAlert("alert2"),
                                             bsAlert("alert3"),
                                             bsAlert("alert4"),
                                             bsAlert("alert5"),
                                             bsAlert("alert6"),
                                             plotOutput("plot1"),
                                             uiOutput(outputId = "download_plot1"),
                                             br(),
                                             br()
                                    ),
                                    tabPanel("Lncrna & Enhancers",
                                             plotOutput("plot2"),
                                             uiOutput(outputId = "download_plot2"),
                                             br(),
                                             plotOutput("plot3"),
                                             uiOutput(outputId = "download_plot3"),
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
                                                      bsButton(inputId = "addRNASeq", label = "Add RNA-Seq analysis", disabled = TRUE),
                                                      br(), br(), bsAlert("addrnaseq_msg_i"), br()),
                                               column(width = 6,
                                                      bsButton(inputId = "runRNASeq",label = "Send analysis request", style = "btn btn-primary", disabled = TRUE),
                                                      br(), br(), bsAlert("runrnaseq_msg_i")))
                                    ),
                                    tabPanel("ChIP-Seq",
                                             h3("ChIP-Seq"),
                                             helpText("We are working with the following datasets : ",
                                                      tags$ul(list(
                                                        tags$li("Targets for MCF7 : CTCF, H3K27ac, H3K36me3, H3K4me3, H3K27me3"),
                                                        tags$li("Targets for K562 : H3K36me3, H3K4me2, CTCF, EZH2"),
                                                        tags$li("Targets for HMEC : H3K4me3, H2AFZ, H3K9me3, H3K36me3, H3K4me3, 
                                             H3K27ac, H3K4me1, H3K79me2, CTCF, H3K9ac, H3K4me2, 
                                             H4K20me1, H3K27me3, EZH2")
                                                      ))),
                                             br(),
                                             br(),
                                             checkboxGroupInput(inputId = "chipseq_analysis",label = "Choose cell types", 
                                                         choices = list("MCF7","K562","HMEC"), selected = "MCF7"),
                                             br(),
                                             br(),
                                             bsButton(inputId = "runCHIPSeq",label = "Send analysis request", style = "btn btn-primary",disabled = TRUE),
                                             br(), br(), bsAlert("runchipseq_msg_i")
                                    )
                                    
                         )
                       )
      )
    )
  ),footer()
))
