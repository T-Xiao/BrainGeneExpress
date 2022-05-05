tabPanel("Multivariate Analysis",
                   fluidRow(
                     #  textOutput('text2'),
                     tags$head(tags$style("#text2{color: black;
                                          font-size: 24px;
                                          font-style: Bold;
                                          }")
                     ),
                     textOutput('text3'),
                     tags$head(tags$style("#text3{color: black;
                                     font-size: 22px;
                                     font-style: Bold;
                                     }")
                     )),
                   sidebarPanel(
                     selectInput("Database","Choose your database",c("TCGA GBM U133a","TCGA GBM RNA-Seq","TCGA GBM Agilent"),width="100%"),
                     selectInput("Target","Choose target gene type",c("Continuous","Binary"),width="600px"),
                     
                     conditionalPanel(condition="input.Database=='TCGA GBM U133a'",
                                      selectInput("Multigene", "Target Gene:", U133aoverallgenenames),
                                      selectInput("Othergenes","Choose other genes",U133aoverallgenenames,selected=NULL,multiple=TRUE,width="600px")
                     ),
                     conditionalPanel(condition="input.Database=='TCGA GBM RNA-Seq'",
                                      selectInput("Multigene2", "Target Gene:", RNAoverallgenenames),
                                      #selectInput("Target","Choose target gene type",c("Continuous","Binary"),width="600px"),
                                      selectInput("Othergenes2","Choose other genes",RNAoverallgenenames,selected=NULL,multiple=TRUE,width="600px")),
                     conditionalPanel(condition="input.Database=='TCGA GBM Agilent'",
                                      selectInput("Multigene3", "Target Gene:", agilentoverallgenenames),
                                      #selectInput("Target","Choose target gene type",c("Continuous","Binary"),width="600px"),
                                      selectInput("Othergenes3","Choose other genes",agilentoverallgenenames,selected=NULL,multiple=TRUE,width="600px")),
                     
                     #selectInput("Endpoint","Choose your endpoint",c("Overall survival","Time-to-recurrence"),width="100%"),
                     selectInput("RCSModel","Choose your model",c("With Restricted Cubic Spline","Without Restricted Cubic Spline"),width="100%"),
                     selectInput("ClinicalFactor","Choose clinical Factors",c("Age","KPS","Gender","Race","G_Cimp_Status","GeneExp_Subtype","Gender","Initial_Pathologic","Radiation_Therapy","Chemo_Therapy"),selected=NULL,multiple=TRUE,width="100%"),
                     
                     selectInput("MissingDataImputation","Missing Data Imputation",c("No imputation","Multiple Imputation"),width="600px"),
                     selectInput("VariableSelection","Choose variable selection methods",c("No selection","Stepwise Selection"),width="600px")
                     ,
                     actionButton("coxcalc", "Run Cox Regression Model"),
                     width=3
                   ), 
                   mainPanel(
                     #numericInput('n', 'Number of obs', 100),
                     tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: relative;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 500%;
               color: #000000;
               background-color: #EBDEF0;
               z-index: 105;
             }
          ")),conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Calculation in Progress...",id="loadmessage")),
                     
                     tabsetPanel(type = "tabs",
                                 tabPanel("Hazard Ratio Table", 
                                          
                                          tableOutput("TCGAtable")),
                                 tabPanel("Hazard Ratio Plot", plotOutput("TCGAplot")),
                                 tabPanel("One and two years Prediction", plotOutput("Predplot1")),
                                 tabPanel('Point Prediction',textInput("inText", "Input text"),tags$img(src="Presentation1.jpg",width = "800px", height = "400px"),plotOutput("Predplot2"))
                     )))