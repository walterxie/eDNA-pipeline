# ui.R

library(shiny)

shinyUI(navbarPage(
  "eDNA pipeline:",
  tabPanel("Upload CM",
           sidebarLayout(
             sidebarPanel(
               fileInput(
                 'cmTaxaFile', 'Upload Community Matrix File : ',
                 accept = c('text/csv',
                            'text/comma-separated-values,text/plain',
                            '.csv')
               ),
               tags$hr(),
               checkboxInput('header', 'Header', TRUE),
               radioButtons('sep', 'Separator',
                            c(
                              Comma = ',',
                              Semicolon = ';',
                              Tab = '\t'
                            ),
                            ','),
               radioButtons(
                 'quote', 'Quote',
                 c(
                   None = '',
                   'Double Quote' = '"',
                   'Single Quote' = "'"
                 ),
                 '"'
               )
             ),
             mainPanel(DT::dataTableOutput('cmTaxaTable'))
           )),
  tabPanel("Taxanomy",
           fluidPage(
             fluidRow(
               column(
                 4,
                 helpText(
                   "Move all taxa to \"Others\" category whose percentage of the total reads <= a given threshold."
                 ),
                 sliderInput(
                   "percThr", "Percentage threshold:", min = 0,
                   max = 1, value = 0.1, step = 0.1, post = "%"
                 ),
                 selectInput("bar_vt", "Values in bar represent:", 
                    choices = c("Species abundance", "Species richness", "OTU abundance", "individuals", "reads", "OTUs"))
               ),
               column(
                 4,
                 #helpText("Only png avaible for web UI. Download pdf for a high resolution image."),
                 numericInput(
                   "legend_nrow", "Set legend row:", 1, min=1
                 ),
                 sliderInput(
                   "zoom", "Zoom in/out", min = 50,
                   max = 200, value = 100, step = 50, post = "%"
                 ),
                 downloadButton('downloadTA', 'Download PDF')
               )
               #column(9, verbatimTextOutput("ta_info"))
             ),
             hr(),
             mainPanel(imageOutput("imageTA"))
           )),
  tabPanel("Summary",
           verbatimTextOutput("summary")),
  navbarMenu(
    "More",
    tabPanel("RCode",
             includeMarkdown("../R/README.md")),
    tabPanel("ReadMe",
             includeMarkdown("README.md")),
    tabPanel("About",
             includeMarkdown("about.md"))
  )
))
