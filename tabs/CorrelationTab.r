tabPanel("Correlation Analysis", 
         sidebarPanel(
                     selectInput("database_corr",
                                 "Select a database",
                                 c("TCGA GBM U133a","TCGA GBM RNA-Seq","TCGA GBM Agilent","TCGA LGG RNA-Seq","IvyGAP"),
                                 selected = NULL,
                                 multiple = F),
                     conditionalPanel(condition="input.Corr_tabs=='One vs Many'",
                                      uiOutput('onemanyone'),
                                      uiOutput('onemanymany')),
                     conditionalPanel(condition="input.Corr_tabs=='One vs All'",
                                      uiOutput('oneallone')),
                                      #select algorithms
                                      # selectInput("Clustering_method","Select the clustering algorithm",
                                      #             c("Hierarchical","DIvisive ANAlysis ",
                                      #               "K-Means","Partition Around Medoids",
                                      #               "Affinity Propagation", "Gaussian Mixture Mode",
                                      #               "Fuzzy C-Means Clustering","HDBSCAN"),multiple = TRUE),
                                      # #conditionalPanel(condition="input.Clustering_method=='Hierachical'",
                                      # selectInput("hc_metric","Select dissmilarity metric",c("euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski")),
                                      # selectInput("hc_method","Select agglomeration method",c("ward.D", "ward.D2", "single", "complete", "average","mcquitty" , "median", "centroid" )),
                                      # selectInput("nk_range","Select how many clusters",c("2 to 7","2 to 8","2 to 9","2 to 10")),
                                      # selectInput("Coef_or_number","Select the method for cutting off",c("Coefficient","Number of Genes Ranking by Absolute Value of Coefficient")),
                                      # conditionalPanel(condition="input.Coef_or_number=='Coefficient'",
                                      #                 sliderInput("cut_off_coef","Please select the cut-off coefficient:", min=-1,max=1,value=c(-1,1),step=0.1)),
                                      # conditionalPanel(condition="input.Coef_or_number=='Number of Genes Ranking by Absolute Value of Coefficient'",
                                      #                                     # sliderInput("cut_off_number","Please select the cut-off number of genes:",min = 0, max = 500,
                                      #                                     #             value = c(0,100))),
                                      #  numericInput("n", "N:", min = 0, max = 500, value = 50))),
                     conditionalPanel(condition="input.Corr_tabs=='Two Gene Sets'",
                                      uiOutput('twosetone'),
                                      uiOutput('twosettwo'),
                                      # selectInput("canodim","Please select the dimensions with a maximum of 3",
                                      #             seq(1,5,1),
                                      #             selected=NULL,multiple=TRUE)
                                      uiOutput('canondimen')
                     ),                                      
                     actionButton("correl", "Run Correlation Analysis")),

         mainPanel(
           
           h3(textOutput(outputId = "caption")),
           
           tabsetPanel(id = "Corr_tabs", 
                       tabPanel(title = "One vs Many", 
                                fluidRow(column(5,plotOutput("ffplot", width =
                                                               "300%", height = "600px")))),
                       tabPanel(title = "One vs All", fluidRow(column(width=12,plotOutput("clusterplot",width="100%",height="600px"))
                       )),
                       
                       tabPanel(title="Two Gene Sets",fluidRow(
                         column(width=12, plotOutput("ccplot", width ="100%", height = "300px"))),
                         fluidRow(
                                  column(width=12,dataTableOutput("cancorrtable")),
                                  column(width=12,dataTableOutput("cancoefx")),
                                  column(width=12,dataTableOutput("cancoefy")),
                                  column(width=12,dataTableOutput("standexvar"))
                         #,
                         #tableOutput("cctable2")
                                                                      ))
           )))

