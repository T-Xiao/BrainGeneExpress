tabPanel(title = "Survival Tree Analysis",
         sidebarPanel(  
           selectInput("Databases",
                       "Choose your database",
                       c("TCGA GBM U133a","TCGA GBM RNA-Seq","TCGA GBM Agilent"),width="600px"),
           conditionalPanel(condition="input.Databases=='TCGA GBM U133a'",
                            selectInput("Unigene", "Target Gene:", U133aoverallgenenames)
           ),
           conditionalPanel(condition="input.Databases=='TCGA GBM RNA-Seq'",
                            selectInput("Unigene2", "Target Gene:", RNAoverallgenenames)),
                            #selectInput("Target","Choose target gene type",c("Continuous","Binary"),width="600px"),
           conditionalPanel(condition="input.Databases=='TCGA GBM Agilent'",
                            selectInput("Unigene3", "Target Gene:", agilentoverallgenenames)),
                            #selectInput("Target","Choose target gene type",c("Continuous","Binary"),width="600px"),
           actionButton("Uni_cox", "Run Univariate Analysis"),
           width=3),
         
         mainPanel(fluidRow(plotOutput("recurplot"),
                            plotOutput("uniplot")
         )))         



  