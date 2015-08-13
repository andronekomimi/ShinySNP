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
    tags$link(rel = "icon", type = "image/x-icon", href = "logo.jpg"),
    tags$script('Shiny.addCustomMessageHandler("myCallbackHandler",
                       function(typeMessage) {console.log(typeMessage)
                           if(typeMessage == 1){
                              $("a:contains(3D Conformation and Regulation)").click();
                           }
                        });')
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
                         hr(),
                         bsCollapse(id = "param", open = c('Advanced parameters'), multiple = FALSE,
                                    bsCollapsePanel('Advanced parameters',uiOutput("highlight"), 
                                                    style = "default")
                         ),
                         bsAlert("addhg_msg_i"),
                         hr(),
                         h3("Add variants"),
                         helpText("Add your variants manually or load a file conaining the variants list"),
                         helpText("File format : ",
                                  tags$ul(
                                    tags$li("Tab-delimited, with 4 named column 'id,start,end,metadata'"),
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
                         h3("Choose the 3D Conformation dataset"),
                         checkboxGroupInput(inputId = "my.dataset",label = NULL, inline = TRUE,
                                            choices = list("local" = "interne", 
                                                           "4DGenome" = "externe"),
                                            selected = c("interne","externe")),
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
                                    tabPanel("3D Conformation and Regulation",
                                             bsAlert("chr_i"),
                                             bsAlert("alert0"),
                                             bsAlert("alert1"),
                                             bsAlert("alert2"),
                                             bsAlert("alert3"),
                                             bsAlert("alert4"),
                                             bsAlert("alert5"),
                                             bsAlert("alert6"),
                                             h4("3D Conformation : Local Dataset"),
                                             plotOutput("plot1"),
                                             uiOutput(outputId = "download_plot1"),
                                             br(),
                                             hr(),
                                             h4("3D Conformation : 4DGenome Dataset"),
                                             plotOutput("plot1b"),
                                             uiOutput(outputId = "download_plot1b"),
                                             br(),
                                             hr(),
                                             h4("Regulation"),
                                             plotOutput("plot2"),
                                             uiOutput(outputId = "download_plot2"),
                                             br(),
                                             h5("LNCRNA Expression levels"),
                                             plotOutput("plot3"),
                                             uiOutput(outputId = "download_plot3"),
                                             br(),
                                             br()
                                    ),
                                    tabPanel("Non-Instant Analysis : RNA-Seq and ChIP-Seq",
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
                                                                   choices = list("Merge experiments" = TRUE,
                                                                                  "One track per experiment" = FALSE
                                                                   ),
                                                                   selected = TRUE)),
                                               column(width = 6,
                                                      bsButton(inputId = "runRNASeq",label = "Run RNA-Seq Analysis", style = "btn btn-primary", disabled = FALSE)
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 12,
                                                      bsAlert("runrnaseq_msg_i"))
                                             ),
                                             br(),
                                             plotOutput("plot4"),
                                             uiOutput(outputId = "download_plot4"),
                                             br(),br(),
                                             hr(),
                                             h3("ChIP-Seq"),
                                             helpText("We are working with the following datasets : ",
                                                      tags$ul(list(
                                                        tags$li("Targets for MCF7 : CTCF, H3K27ac, H3K36me3, H3K4me3, H3K27me3"),
                                                        tags$li("Targets for K562 : H3K36me3, H3K4me2, CTCF, EZH2"),
                                                        tags$li("Targets for HMEC : H2AFZ, H3K9me3, H3K36me3, H3K4me3, 
                                             H3K27ac, H3K4me1, H3K79me2, CTCF, H3K9ac, H3K4me2, 
                                             H4K20me1, H3K27me3, EZH2")
                                                      ))),
                                             br(),
                                             fluidRow(
                                               column(width = 6,
                                                      radioButtons(inputId = "chipseq_cell",label = "Choose cell type", 
                                                                   choices = list("MCF7" = "mcf7" ,"K562" = "k562","HMEC" = "hmec"), selected = "mcf7")),
                                               column(width = 6,
                                                      bsButton(inputId = "runCHIPSeq",label = "Run ChIP-Seq Analysis", style = "btn btn-primary", disabled = TRUE)
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 12,
                                                      bsAlert("runchipseq_msg_i"))
                                             ),
                                             br(),
                                             plotOutput("plot5"),
                                             uiOutput(outputId = "download_plot5"),
                                             br()
                                    )
                                    
                         )
                       )
      )
    )
  ),footer()
))
