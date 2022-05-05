shinyServer(function(input, output, session) {
  ### Correlation tab:
  # Load correlation databases according to users' choice:
  correlation_database <- eventReactive(input$correl,
                                        isolate({
                                          if (input$database_corr == "TCGA GBM U133a") {
                                            corr_database = TCGA_GBM_u133a_overall[, 1:12042]
                                          } else if (input$database_corr == "TCGA GBM RNA-Seq") {
                                            corr_database = TCGA_GBM_rna_seq_overall[, 1:20028]
                                          } else if (input$database_corr == "TCGA GBM Agilent") {
                                            corr_database = TCGA_GBM_agilent_overall[, 1:17814]
                                          } else if (input$database_corr == "TCGA LGG RNA-Seq") {
                                            corr_database = TCGA_LGG_overall[, 1:20225]
                                          } else if (input$database_corr == "IvyGAP") {
                                            corr_database = Ivygap[, 1:25873]
                                          }
                                          return(corr_database)
                                        }))
  # Load different gene name lists for different data sets:
  datasetInput <- reactive({
    switch(
      input$database_corr,
      "TCGA GBM U133a" = U133aoverallgenenames,
      "TCGA GBM RNA-Seq" = RNAoverallgenenames,
      "TCGA GBM Agilent" = agilentoverallgenenames,
      "TCGA LGG RNA-Seq" = lggoverallgenenames,
      "IvyGAP" = ivygapgenenames
    )
  })
  ## One vs many correlation dataset input:
  output$onemanyone = renderUI({
    mydata = datasetInput()
    selectInput('onemanyone', 'Select One Gene', mydata)
  })
  output$onemanymany = renderUI({
    mydata = datasetInput()
    selectInput(
      'onemanymany',
      'Select Several Genes',
      mydata,
      selected = NULL,
      multiple = TRUE
    )
  })
  ## one vs all correlation dataset input:
  output$oneallone = renderUI({
    mydata = datasetInput()
    #adding a div to 
    div(selectInput("oneall", "Please select one genes:", mydata),
    selectInput("type_of_coef","Please select the type of correlation coefficient",c("pearson","spearman","kendall")),
   # sliderInput("cut_off_coef","Please select the cut-off correlation coefficient:", min=-1,max=1,value=c(-1,1),step=0.1),
    selectInput("method_of_ranking","Please select the ranking method of coefficient",c("Positive","Negative","Mixed")),
    numericInput("n_genes", "Enter Number of Genes (1 to 500)", min = 1, max = 500, value = 50),
    #select algorithms
    selectInput("Clustering_method","Select the clustering algorithm(s)",
                c("Hierarchical","DIvisive-ANAlysis",
                  "K-Means","Partition-Around-Medoids",
                  "Affinity-Propagation", "Gaussian-Mixture-Mode",
                  "Fuzzy-C-Means-Clustering","HDBSCAN"),multiple = TRUE,selected="Hierarchical"),
    selectInput("hc_metric","Select dissmilarity metric",c("euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski")),
    selectInput("hc_method","Select agglomeration method",c("ward.D", "ward.D2", "single", "complete", "average","mcquitty" , "median", "centroid" )),
    selectInput("nk_range","Select range of clusters",c("2 to 7","2 to 8","2 to 9","2 to 10")),
    selectInput("n_best","Select the top n algorithms to keep after trimming off the poor performing ones using Rank Aggregation",c(rep(1:8,1))))
    #                                     # sliderInput("cut_off_number","Please select the cut-off number of genes:",min = 0, max = 500,
    #                                     #             value = c(0,100))),
  })
  ## two sets correlation dataset input:
  output$twosetone = renderUI({
    mydata = datasetInput()
    selectInput(
      "Two_sets_one",
      "Please select the 1st list of genes:",
      mydata,
      selected = NULL,
      multiple = TRUE
    )
  })
  output$twosettwo = renderUI({
    mydata = datasetInput()
    selectInput(
      "Two_sets_two",
      "Please select the 2nd list of genes:",
      mydata,
      selected = NULL,
      multiple = TRUE
    )
  })
  output$canondimen = renderUI({
    selectizeInput(
      "canodim",
      paste0("Please select up to 3 dimensions:"),
      c(seq(1, min(
        length(req(input$Two_sets_one)), length(req(input$Two_sets_two))
      ), by = 1)),
      selected = 1,
      multiple = TRUE,
      options = list(maxItems = 3)
    )
  })
  ## one vs many plot:
  output$ffplot <- renderPlot({
    corr_database = correlation_database()
    one_gene = input$onemanyone
    many_gene = strsplit(as.character(input$onemanymany), " ")
    ngenelist = length(unlist(many_gene))
    #correlation coefficient:
    # corr_eqn <- function(x, y, digits = 2) {
    #   corr_coef <- round(cor(x, y), digits = digits)
    #   paste("r", "=", corr_coef)
    # }
    # p_eqn <- function(x, y, digits = 4) {
    #   p <- cor.test(x, y, method = "pearson")$p.value
    #   txt2 <-
    #     format(c(p, 0.123456789), digits = digits)[1]
    #   paste("p= ", txt2, sep = "")
    # }
    plotgene = vector("list", length = ngenelist)
    #  print (plotgene)
    for (i in 1:ngenelist) {
      #local({
      gene_name = unlist(many_gene)[i]
      print(gene_name)
      pp <<-
        ggplot(corr_database,
               aes_string(corr_database[, one_gene], corr_database[, gene_name])) +
        geom_point(shape = 1) +
        geom_smooth(method = lm) +
        ylab(one_gene) +
        xlab(gene_name) +
        # ylim(0,1.1e5) +
        geom_cor(method = "pearson", ypos = 1e5) 
        # geom_text(
        #   y = max(corr_database[, one_gene]),
        #   x = max(corr_database[, gene_name]) - 0.5,
        #   label = corr_eqn(corr_database[, gene_name], corr_database[, one_gene])
        # ) +
        # geom_text(
        #   y = max(corr_database[, one_gene]) - 0.5,
        #   x = max(corr_database[, gene_name]) - 0.5,
        #   label = p_eqn(corr_database[, gene_name], corr_database[, one_gene])
        # )
      plotgene[[i]] = pp
      #print(i)
      #print(g)
      # assign(unlist(many_gene)[i],pp,pos=1)})
    }
    # print (length(plotgene))
    plotgene
    grid.arrange(grobs = plotgene, ncol = 3)
  })
  ##one vs all plots:
  one_vs_all <- eventReactive(input$correl,
                              isolate({
                                corr_database = correlation_database()
                                one_all_gene = input$oneall
                                correlation = data.frame(cor(corr_database[,one_all_gene], corr_database,method=input$type_of_coef))
                                correlation_tr = data.frame(t(correlation))
                                #depending on "positive", "negative" and "all",different correlation order dataframe:
                                if (input$method_of_ranking=="Positive"){
                                  correlation_tr_order = data.frame(correlation_tr[order(-(correlation_tr)), , drop = FALSE])
                                } else if (input$method_of_ranking=="Negative"){
                                  correlation_tr_order = data.frame(correlation_tr[order(correlation_tr), , drop = FALSE])
                                } else if (input$method_of_ranking=="Mixed"){
                                correlation_tr_order = data.frame(correlation_tr[order(-abs(correlation_tr)), , drop = FALSE])
                                }
                               
                                if (input$n_genes > 500){
                                  showNotification("Number of genes is set to 500",type="warning")
                                  ngenes <- 500
                                } else if (input$n_genes < 1){
                                  showNotification("Number of genes is set to 1",type="warning")
                                  ngenes <- 1
                                  }else {ngenes <- input$n_genes}
                                
                                heatmapdf = as.data.frame(t(corr_database[, rownames(correlation_tr_order)[1:ngenes]]))
                                return(heatmapdf)
                              }))
  output$heatmapdf <- renderTable({heatmap1 = one_vs_all()
                                   return (heatmap1)})

  output$clusterplot <- renderPlot({
    heatmap1 = one_vs_all()
    # convert strings to readable range of nk:
    if (input$nk_range=="2 to 7")
    {nkrange=2:7
      }else if (input$nk_range=="2 to 8" )
      {nkrange=2:8
      }else if (input$nk_range=="2 to 9")
      {nkrange=2:9
      }else {nkrange=2:10}
    # transform the algorithms to what dice can read:
    algorithms = strsplit(as.character(input$Clustering_method), " ")
    algorithms_replaced<-replace(algorithms,c(which(algorithms%in%"Hierarchical"),which(algorithms%in%"K-Means"),
                             which(algorithms%in%"DIvisive-ANAlysis"),which(algorithms%in%"Partition-Around-Medoids"),
                             which(algorithms%in%"Affinity-Propagation"),which(algorithms%in%"Gaussian-Mixture-Mode"),
                             which(algorithms%in%"Fuzzy-C-Means-Clustering"),which(algorithms%in%"HDBSCAN")),
                             c("hc","km","diana","pam","ap","gmm","cmeans","hdbscan"))
    
    diceR_result <- dice(as.matrix(heatmap1),nk=nkrange,algorithms =unlist(algorithms_replaced), reps=10,
                                        #hc.method = input$hc_method,
                                         #distance = input$hc_metric, 
                                         cons.funs = c( "majority", "CSPA","LCE"), 
                                         sim.mat=c("cts"),
                                         prep.data = c("full"),
                                          min.var = 1, seed = 1,
                                         trim = TRUE, reweigh = TRUE, evaluate = FALSE,n=6,
                                         plot = FALSE, ref.cl = NULL, progress = TRUE)
    clus_gene_exp <- merge(heatmap1,diceR_result$clusters[,1],by="row.names")
    clus_gene_exp_final <- clus_gene_exp %>% .[order(.$y),] %>% set_rownames(.$Row.names) %>% dplyr::select(-Row.names)
    pheatmap(clus_gene_exp_final,cluster_rows = FALSE,cluster_cols = FALSE,scale="column",show_colnames=FALSE,fontsize_row=4)
    
  })
  ## Two gene sets canonical correlation:
  ## We need both cancor file and
  #Create a cancor file for table and plot; a CanCorr.r file for table
  #Extract Two gene sets:
  twogenesets <- reactive({
    TCGA = correlation_database()
    split3 = strsplit(input$Two_sets_one, " ")
    #  create Y from genelist 3 for canonical correlation coefficient:
    split4 = strsplit(input$Two_sets_two, " ")
    return (list(split3, split4))
  })
  # Generate candisc::cancor file, yacca::cca file and yacca:F.test.cca file to
  # be used:
  cancorreact <- reactive({
    corr_database = correlation_database()
    X <- twogenesets()[[1]]
    Y <- twogenesets()[[2]]
    candisc_cca <-
      cancor(corr_database[unlist(X)],
             corr_database[unlist(Y)],
             set.names = c("Gene Set 1", "Gene Set 2"))
    yacca_cancor <-
      cca(corr_database[unlist(X)],
          corr_database[unlist(Y)],
          xscale = TRUE,
          yscale = TRUE)
    #Wilkscca<-Wilks(cca)
    yacca_Ftest <- F.test.cca(yacca_cancor)
    return (list(candisc_cca, yacca_cancor, yacca_Ftest))
  })
  #Render (1) a datatable with standardized canonical correlation, p.value, etc.
  output$cancorrtable <- renderDataTable({
    x <- cancorreact()[[3]]
    signif <- c(ifelse(
      x$p.value < 0.0001,
      "****",
      ifelse(
        x$p.value < 0.001,
        "***",
        ifelse(x$p.value < 0.01, "**",
               ifelse(x$p.value < 0.05, "*", "ns"))
      )
    ))
    datatable(
      data.frame(
        Corr = x$corr,
        CanRSq = x$corr ^ 2,
        Statistic = x$statistic,
        parameter = x$parameter,
        p.value = x$p.value,
        signficance = signif
      ) %>%
        set_colnames(
          c(
            "Stand.Can.Corr.",
            "Squared Can.Corr.",
            "F",
            "Num df",
            "Den df",
            "Pr(>F)",
            "Sig. Level"
          )
        ),
      options = list(pageLength = 10, lengthChange = FALSE),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: center;',
        htmltools::hr('Table 1: Standardized Canonical Correlation Table')
      )
    ) %>%
      formatRound(c(1:7), 3) %>%
      formatStyle(columns = c(1:7), 'text-align' = 'center')
  })
  #Render (2) a datatable of x canonical coefficient with yacca::cca
  output$cancoefx <- renderDataTable({
    x <- cancorreact()[[2]]
    datatable(
      df <-
        data.frame(x$xcoef),
      options = list(pagelength = 10, lengthChange = FALSE),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: center;',
        htmltools::hr('Table 2:Gene Set 1 Canonical Coefficients')
      )
    ) %>%
      formatRound(c(1:ncol(df)), 3) %>%
      formatStyle(columns = c(1:ncol(df)), 'text-align' = 'right')
  })
  #Render (3) a datatable of y canonical coefficient with yacca::cca
  output$cancoefy <- renderDataTable({
    x <- cancorreact()[[2]]
    datatable(
      df <-
        data.frame(x$ycoef),
      options = list(pagelength = 10, lengthChange = FALSE),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: center;',
        htmltools::hr('Table 3: Gene Set 2 Canonical Coefficients')
      )
    ) %>%
      formatRound(c(1:ncol(df)), 3) %>%
      formatStyle(columns = c(1:ncol(df)), 'text-align' = 'right')
  })
  #Render (4) a datatable for standardized explained variance:
  output$standexvar <- renderDataTable({
    x <- cancorreact()[[2]]
    sketch = htmltools::withTags(table(class = 'display',
                                       # style = 'text-align: center;',
                                       thead(tr(
                                         th(
                                           style = 'text-align: center;',
                                           rowspan = 2,
                                           'Canonical Variates'
                                         ),
                                         th(style = 'text-align: center;', colspan = 2, 'Gene Set 1 '),
                                         th(style = 'text-align: center;', colspan = 2, 'Gene Set 2')
                                       ),
                                       tr(
                                         lapply(rep(
                                           c('Proportion (Own)', 'Proportion (Opposite)'), 2
                                         ), th,
                                         style = 'text-align: center;')
                                       ))))
    # x<-cancorreact()[[2]]
    # vec<-c()
    # sumxvrd<-0
    # vec2<-c()
    # sumxcanvad<-0
    # for (i in 1:length(x$corr)){
    #   sumxvrd<-x$xvrd[i]+sumxvrd
    #   vec=c(vec,sumxvrd)
    # }
    # for (i in 1:length(x$corr)){
    #   sumxcanvad<-x$xcanvad[i]+sumxcanvad
    #   vec2=c(vec2,sumxcanvad)
    # }
    k <- length(x$corr)
    
    datatable(
      df <- data.frame(
        variates = paste("CV", 1:k),
        own_can_va1 = x$xcanvad,
        #sumcanvad=vec2,
        opp_can_va1 = x$xvrd,
        own_can_va2 = x$ycanvad,
        opp_can_va2 = x$yvrd
      ),
      container = sketch,
      rownames = FALSE
    ) %>% formatRound(c(2:ncol(df)), 3) %>%
      formatStyle(columns = c(1:ncol(df)), 'text-align' = 'center')
    
  })
  #Render (5) cca plot
  output$ccplot <- renderPlot({
    cca <- cancorreact()[[1]]
    par(mfrow = c(1, 3))
    if (length(input$canodim) == 3) {
      plot(cca, which = strtoi(input$canodim[1]))
      plot(cca, which = strtoi(input$canodim[2]))
      plot(cca, which = strtoi(input$canodim[3]))
    }
    else if (length(input$canodim) == 2) {
      plot(cca, which = strtoi(input$canodim[1]))
      plot(cca, which = strtoi(input$canodim[2]))
    }
    else if (length(input$canodim) == 1) {
      plot(cca, which = strtoi(input$canodim[1]))
    }
  })
  ### Multivariable Tab
  calcultatecox <-eventReactive(input$coxcalc,
                                isolate({if (input$Database == "TCGA GBM RNA-Seq") {
                                  TCGA = TCGA_GBM_rna_seq_overall
                                  target = TCGA[, input$Multigene2]
                                } else if (input$Database == "TCGA GBM U133a") {
                                  TCGA = TCGA_GBM_u133a_overall
                                  target = TCGA[, input$Multigene]
                                } else if (input$Database == "TCGA GBM Agilent") {
                                  TCGA = TCGA_GBM_agilent_overall
                                  target = TCGA[, input$Multigene3]}
                                  #print (TCGA[,1])
                                  return(list(TCGA, target))
                                }))
  output$text3 <-
    renderText("Cox Proportional Hazards Model")
  sfreact <- reactive({
    input$coxcalc
    isolate({TCGA = calcultatecox()[[1]]
    #print (TCGA[,1])
    TCGA$target_gene = calcultatecox()[[2]]
    if (input$Database == "TCGA GBM RNA-Seq") {
      targetname = input$Multigene2
      Othergenes = input$Othergenes2
      TCGAname = rep("TCGA", 20028)
      medianSelected = median(TCGA$target_gene)
      
      for (i in 1:172) {
        if (TCGA$target_gene[i] < medianSelected) {
          TCGA$Target_Gene_Stratified[i] = "LOW"
        }
        else {
          TCGA$Target_Gene_Stratified[i] = "HIGH"
        }
      }
    } else if (input$Database == "TCGA GBM U133a") {
      targetname = input$Multigene
      Othergenes = input$Othergenes
      TCGAname = rep("TCGA", 12042)
      
      medianSelected = median(TCGA$target_gene)
      for (i in 1:539) {
        if (TCGA$target_gene[i] < medianSelected) {
          TCGA$Target_Gene_Stratified[i] = "LOW"
        }
        else {
          TCGA$Target_Gene_Stratified[i] = "HIGH"
        }
      }
    }  else if (input$Database == "TCGA GBM Agilent") {
      targetname = input$Multigene3
      Othergenes = input$Othergenes3
      TCGAname = rep("TCGA", 17814)
      
      medianSelected = median(TCGA$target_gene)
      for (i in 1:585) {
        if (TCGA$target_gene[i] < medianSelected) {
          TCGA$Target_Gene_Stratified[i] = "LOW"
        }
        else {
          TCGA$Target_Gene_Stratified[i] = "HIGH"
        }
      }
    }
    print(targetname)
    #sum(TCGA_GBM_rna_seq_overall$RNF17)
    
    split_factor = strsplit(as.character(input$ClinicalFactor), " ")
    split_cont = strsplit(as.character(Othergenes), " ")
    
    if (is.element('Age', split_factor)) {
      split_factor = split_factor[split_factor != "Age"]
      split_cont = list.append(split_cont, 'Age')
    } else {
      split_factor = split_factor
      split_cont = split_cont
    }
    # duplicate target gene: delete one
    if (input$Target == "Continuous") {
      split_cont = unique(list.append(split_cont, targetname))
    }
    else if (input$Target == "Binary") {
      split_factor = list.append(split_factor, "Target_Gene_Stratified")
    }
    dd <<- datadist(TCGA)
    options(datadist = 'dd')
    if (input$RCSModel == "With Restricted Cubic Spline") {
      if (input$MissingDataImputation == "Multiple Imputation") {
        #combine to do imputation:
        impute_for_all_list = list.append(split_cont, split_factor)
        form_impute <-
          paste("~I(Time_To_Event)+ Event+",
                paste(unlist(impute_for_all_list), collapse = "+"))
        set.seed(133344)
        ## test if there will be an error:
        if (is.error(
          try(aregImpute(
            as.formula(form_impute),
            data = TCGA,
            nk =3,
            n.impute = 10
          )))==TRUE){a <-
            aregImpute(
              as.formula(form_impute),
              data = TCGA,
              nk =0,
              n.impute = 10)
        } else {a <-
          aregImpute(
            as.formula(form_impute),
            data = TCGA,
            nk =3,
            n.impute = 10)}
        # 
        # tryCatch(
        #   a <-
        #       aregImpute(
        #         as.formula(form_impute),
        #         data = TCGA,
        #        nk =3,
        #         n.impute = 10
        #       ),
        #   
        #   return(a <-
        #            aregImpute(
        #              as.formula(form_impute),
        #              data = TCGA,
        #              nk =0,
        #              n.impute = 10))
        #   
        # , error = function(warning) {
        #   message("My message is here!")
        #   return(NA)}
        #       )
        #   return (a)
        # }
        # )
        # 
        # a <-
        #   aregImpute(
        #     as.formula(form_impute),
        #     data = TCGA,
        #    nk =3,
        #     n.impute = 10
        #   )
        # print (a)
        # do formula only for factor=0 or both:
        if (length(split_cont) == 0) {
          formu_new <- paste("Surv(Time_To_Event,Event)~",
                             paste(
                               unlist(split_factor),
                               collapse = "+",
                               sep = ""
                             ))
        } else if (length(split_cont) != 0) {
          if (length(split_factor) == 0) {
            formu <-
              paste(
                "Surv(Time_To_Event, Event)~rcs(",
                paste(
                  unlist(split_cont),
                  collapse = ",3)+rcs(",
                  sep = ""
                ),
                ")"
              )
          }
          else {
            formu <-
              paste(
                "Surv(Time_To_Event,Event)~rcs(",
                paste(
                  unlist(split_cont),
                  collapse = ",3)+rcs(",
                  sep = ""
                ),
                ")+",
                paste(
                  unlist(split_factor),
                  collapse = "+",
                  sep = ""
                )
              )
          }
          #use g formula to test if rcs is needed or not:
          g <- fit.mult.impute(as.formula(formu),
                               cph,
                               a,
                               family = cph,
                               data = TCGA)
          anov <- anova(g)
          rcslist <- list()
          nonlist <- list()
          for (i in seq(2, (length(split_cont)) * 2, 2)) {
            pval = anov[i, 3]
            if (pval < 0.05)
            {
              rcslist = c(rcslist, rownames(anov)[i - 1])
            }
            else
            {
              nonlist = c(nonlist, rownames(anov)[i - 1])
            }
          }
          # print (rcslist)
          # print (nonlist)
          ## add
          nonlist = list.append(nonlist, split_factor)
          ## final formula: three situations:
          if (length(nonlist) == 0) {
            formu_new <- paste(
              "Surv(Time_To_Event, Event)~rcs(",
              paste(
                unlist(rcslist),
                collapse = ",3)+rcs(",
                sep = ""
              ),
              ")"
            )
          } else if (length(rcslist) == 0) {
            formu_new <- paste("Surv(Time_To_Event, Event)~",
                               paste(
                                 unlist(nonlist),
                                 collapse = "+",
                                 sep = ""
                               ))
          } else {
            formu_new = paste(
              "Surv(Time_To_Event,Event)~rcs(",
              paste(
                unlist(rcslist),
                collapse = ",3)+rcs(",
                sep = ""
              ),
              ")+",
              paste(
                unlist(nonlist),
                collapse = "+",
                sep = ""
              )
            )
          }
        }
        ## model selection or not:
        #model selection-- stepwise vs no selection
        if (input$VariableSelection == "Stepwise Selection") {
          g <- fit.mult.impute(as.formula(formu_new),
                               cph,
                               a,
                               family = cph,
                               data = TCGA)
          print (g)
          if (input$Target == "Continuous") {
            scopy = as.formula(paste("~", targetname))
          }
          else {
            scopy = as.formula(paste("~", "Target_Gene_Stratified"))
          }
          
          gfit <-
            stepAIC(
              g,
              direction = "both",
              scope = list(lower = scopy),
              k = 2
            )
          
          print (gfit$formula)
          final_formu = gfit$formula
          noofrcs = grep("\\<rcs\\>", deparse(final_formu))
          # rcslist=list()
          # for (i in 0:noofrcs){
          #   rcslist=list.append(gfit$Design$name[i],rcslist)
          # }
        }
        else if (input$VariableSelection == "No selection") {
          final_formu = formu_new
          # rcslist=rcslist
        }
      }
      
      else if (input$MissingDataImputation == "No imputation") {
        if (length(split_cont) == 0) {
          formu_new <- paste("Surv(Time_To_Event,Event)~",
                             paste(
                               unlist(split_factor),
                               collapse = "+",
                               sep = ""
                             ))
        } else if (length(split_cont) != 0) {
          if (length(split_factor) == 0) {
            formu <-
              paste(
                "Surv(Time_To_Event, Event)~rcs(",
                paste(
                  unlist(split_cont),
                  collapse = ",3)+rcs(",
                  sep = ""
                ),
                ")"
              )
          }
          else  {
            formu <-
              paste(
                "Surv(Time_To_Event,Event)~rcs(",
                paste(
                  unlist(split_cont),
                  collapse = ",3)+rcs(",
                  sep = ""
                ),
                ")+",
                paste(
                  unlist(split_factor),
                  collapse = "+",
                  sep = ""
                )
              )
          }
          
          fcph <-
            cph(as.formula(formu),
                data = TCGA,
                surv = TRUE)
          anov <- anova(fcph)
          rcslist <- list()
          nonlist <- list()
          for (i in seq(2, (length(split_cont)) * 2, 2)) {
            pval = anov[i, 3]
            if (pval < 0.05)
            {
              rcslist = c(rcslist, rownames(anov)[i - 1])
            }
            else
            {
              nonlist = c(nonlist, rownames(anov)[i - 1])
            }
          }
          nonlist <- list.append(nonlist, split_factor)
          print (nonlist)
          
          
          
          #print (rcslist)
          #print (length(rcslist[[1]]))
          #nonlist=list.append(nonlist,)
          # 2018/5/8: Stop here
          
          ## final formula: three situations:
          if (length(nonlist[[1]]) == 0) {
            formu_new <- paste(
              "Surv(Time_To_Event, Event)~rcs(",
              paste(
                unlist(rcslist),
                collapse = ",3)+rcs(",
                sep = ""
              ),
              ")"
            )
          }
          else if (length(rcslist) == 0) {
            formu_new <- paste("Surv(Time_To_Event, Event)~",
                               paste(
                                 unlist(nonlist),
                                 collapse = "+",
                                 sep = ""
                               ))
          }
          
          else if (length(rcslist[[1]]) != 0 &
                   length(nonlist[[1]]) != 0) {
            formu_new = paste(
              "Surv(Time_To_Event,Event)~rcs(",
              paste(
                unlist(rcslist),
                collapse = ",3)+rcs(",
                sep = ""
              ),
              ")+",
              paste(
                unlist(nonlist),
                collapse = "+",
                sep = ""
              )
            )
            print (formu_new)
            
          }
        }
        
        ## model selection or not:
        #model selection-- stepwise vs no selection
        if (input$VariableSelection == "Stepwise Selection") {
          fcph_new <-
            cph(as.formula(formu_new),
                data = na.omit(TCGA),
                surv = TRUE)
          
          fcph_fit <-
            stepAIC(
              fcph_new,
              direction = "both",
              scope = list(lower = formu_new),
              k = 2
            )
          ##20180427 stop here
          #
          #       if (is.null(fcph_fit)==TRUE){
          #   final_formu="ERROR"
          # }else {
          final_formu = fcph_fit$sformula
        }
        
        else if (input$VariableSelection == "No selection") {
          final_formu = formu_new
        }
      }
    }
    else if (input$RCSModel == "Without Restricted Cubic Spline") {
      if (input$MissingDataImputation == "Multiple Imputation") {
        #combine to do imputation:
        impute_for_all_list = list.append(split_cont, split_factor)
        form_impute <-
          paste("~I(Time_To_Event)+ Event+",
                paste(unlist(impute_for_all_list), collapse = "+"))
        set.seed(133344)
        a <-
          aregImpute(
            as.formula(form_impute),
            data = TCGA,
            nk = 3,
            n.impute = 10
          )
        formu_new <- paste("Surv(Time_To_Event,Event)~",
                           paste(
                             unlist(impute_for_all_list),
                             collapse = "+",
                             sep = ""
                           ))
        
        ## model selection or not:
        #model selection-- stepwise vs no selection
        if (input$VariableSelection == "Stepwise Selection") {
          g <- fit.mult.impute(as.formula(formu_new),
                               cph,
                               a,
                               family = cph,
                               data = TCGA)
          if (input$Target == "Continuous") {
            scopy = as.formula(paste("~", targetname))
          }
          else {
            scopy = as.formula(paste("~", "Target_Gene_Stratified"))
          }
          
          gfit <-
            stepAIC(
              g,
              direction = "both",
              scope = list(lower = scopy),
              k = 2
            )
          print (gfit$formula)
          
          final_formu = gfit$formula
          
        }
        else if (input$VariableSelection == "No selection") {
          final_formu = formu_new
        }
      }
      else if (input$MissingDataImputation == "No imputation") {
        nonlist = list.append(split_cont, split_factor)
        formu_new <- paste("Surv(Time_To_Event,Event)~",
                           paste(
                             unlist(nonlist),
                             collapse = "+",
                             sep = ""
                           ))
        if (input$VariableSelection == "Stepwise Selection") {
          fcph_new <-
            cph(as.formula(formu_new),
                data = na.omit(TCGA),
                surv = TRUE)
          fcph_fit <-
            stepAIC(
              fcph_new,
              direction = "both",
              scope = list(lower = formu_new),
              k = 2
            )
          final_formu = fcph_fit$sformula
        }
        
        else if (input$VariableSelection == "No selection") {
          final_formu = formu_new
        }
      }
    }
    #assign variables to a list for prediction
    # varlist<-list()
    # for (i in 3:length(all.vars(final_formu))){
    #   varlist<-list.append(varlist,all.vars(final_formu)[i])
    # }
    print (final_formu)
    varlist <-
      list(all.vars(as.formula(final_formu))[3:length(all.vars(as.formula(final_formu)))])
    
    return(list(
      final_formu = final_formu,
      a = a,
      varlist = varlist
    ))
    # else if (input$Target == "continuous"){
    # return (list(final_formu,a))}
    
    })
  })
  output$TCGAtable <- renderTable({
    ## if rcs is selected:
    #if input$Model==With Restricted Cubic Spline
    input$coxcalc
    isolate({
      final_formu = sfreact()$final_formu
      # rcslist = sfreact()[[2]]
      #print (rcslist)
      a = sfreact()$a
      varlist = sfreact()$varlist
      
      #validate( need(input$VariableSelection =="No selection","No variable left after model selection!") )
      
      TCGA = calcultatecox()[[1]]
      #print (TCGA[,1])
      TCGA$target_gene = calcultatecox()[[2]]
      
      # TCGA$target_gene = calcultatecox()
      medianSelected = median(TCGA$target_gene)
      if (input$Database == "TCGA GBM RNA-Seq") {
        targetname = input$Multigene2
        Othergenes = input$Othergenes2
        TCGAname = rep("TCGA", 20028)
        medianSelected = median(TCGA$target_gene)
        
        for (i in 1:172) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      } else if (input$Database == "TCGA GBM U133a") {
        targetname = input$Multigene
        Othergenes = input$Othergenes
        TCGAname = rep("TCGA", 12042)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:539) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }  else if (input$Database == "TCGA GBM Agilent") {
        targetname = input$Multigene3
        Othergenes = input$Othergenes3
        TCGAname = rep("TCGA", 17814)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:585) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }
      # for (i in 1:539) {
      #   if (TCGA$target_gene[i] < medianSelected) {
      #     TCGA$Target_Gene_Stratified[i] = "LOW"
      #   }
      #   else {
      #     TCGA$Target_Gene_Stratified[i] = "HIGH"
      #   }}
      #
      # if (input$Endpoint=="Time-to-recurrence"){
      #   TCGA$Event=TCGA$recurrence
      #   TCGA$Time_To_Event=TCGA$time_to_recurrence
      # }
      # else if (input$Endpoint=="Overall survival"){
      #   TCGA$Event=TCGA$Event
      #   TCGA$Time_To_Event=TCGA$Time_To_Event_TO_Event}
      dd <<- datadist(TCGA)
      options(datadist = 'dd')
      if (input$MissingDataImputation == "Multiple Imputation") {
        sf <- fit.mult.impute(as.formula(final_formu),
                              cph,
                              a,
                              family = cph,
                              data = TCGA)
      }
      else if (input$MissingDataImputation == "No imputation") {
        sf <- cph(as.formula(final_formu),
                  data = TCGA,
                  surv = TRUE)
      }
      # summary of :
      ssf = summary(sf)
      # check how many variables:
      noofvariables = dim(ssf)[1]
      # calculate p value:
      pvalue = c()
      for (i in 1:(noofvariables / 2)) {
        pvalue = round(c(pvalue, exp(
          -0.717 * abs(ssf[2 * i - 1, 4] / ssf[2 * i - 1, 5]) - 0.416 * abs(ssf[2 *
                                                                                  i - 1, 4] / ssf[2 * i - 1, 5]) ^ 2
        )), 4)
      }
      ##parameter
      variablelist = c()
      for (i in seq(1, noofvariables, 2)) {
        variablelist = c(variablelist, rownames(ssf)[i])
      }
      
      #Hazard ratio
      finaldf = as.data.frame(matrix(ncol = 5, nrow = (noofvariables / 2)))
      colnames(finaldf) = c("Parameters",
                            "Hazard Ratio",
                            "Lower 95%",
                            "Upper 95%",
                            "P value")
      for (i in 1:(noofvariables / 2)) {
        finaldf[i, 1] = variablelist[i]
        finaldf[i, 2] = ssf[2 * i, 4]
        finaldf[i, 3] = ssf[2 * i, 6]
        finaldf[i, 4] = ssf[2 * i, 7]
        finaldf[i, 5] = pvalue[i]
      }
      finaldf})
  })
  output$TCGAplot <- renderPlot({ 
    #if input$Model==With Restricted Cubic Spline
    input$coxcalc
    isolate({
      final_formu = sfreact()$final_formu
      
      # rcslist = sfreact()[[2]]
      #print (rcslist)
      a = sfreact()$a
      varlist = sfreact()$varlist
      TCGA = calcultatecox()[[1]]
      #print (TCGA[,1])
      TCGA$target_gene = calcultatecox()[[2]]
      #TCGA$target_gene = calcultatecox()
      if (input$Database == "TCGA GBM RNA-Seq") {
        targetname = input$Multigene2
        Othergenes = input$Othergenes2
        TCGAname = rep("TCGA", 20028)
        medianSelected = median(TCGA$target_gene)
        
        for (i in 1:172) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      } else if (input$Database == "TCGA GBM U133a") {
        targetname = input$Multigene
        Othergenes = input$Othergenes
        TCGAname = rep("TCGA", 12042)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:539) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }  else if (input$Database == "TCGA GBM Agilent") {
        targetname = input$Multigene3
        Othergenes = input$Othergenes3
        TCGAname = rep("TCGA", 17814)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:585) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }
      
      dd <<- datadist(TCGA)
      options(datadist = 'dd')
      if (input$MissingDataImputation == "Multiple Imputation") {
        sf <- fit.mult.impute(as.formula(final_formu),
                              cph,
                              a,
                              family = cph,
                              data = TCGA)
      } else if (input$MissingDataImputation == "No imputation") {
        sf = cph(as.formula(final_formu),
                 data = TCGA,
                 surv = TRUE)
      }
      ssf = summary(sf)
      plot(ssf)
    })
  })
  output$Predplot1 <- renderPlot({
    #TCGA$target_gene = calcultatecox()
    input$coxcalc
    isolate({
      TCGA = calcultatecox()[[1]]
      #print (TCGA[,1])
      TCGA$target_gene = calcultatecox()[[2]]
      medianSelected = median(TCGA$target_gene)
      
      if (input$Database == "TCGA GBM RNA-Seq") {
        targetname = input$Multigene2
        Othergenes = input$Othergenes2
        TCGAname = rep("TCGA", 20028)
        medianSelected = median(TCGA$target_gene)
        
        for (i in 1:172) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      } else if (input$Database == "TCGA GBM U133a") {
        targetname = input$Multigene
        Othergenes = input$Othergenes
        TCGAname = rep("TCGA", 12042)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:539) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }  else if (input$Database == "TCGA GBM Agilent") {
        targetname = input$Multigene3
        Othergenes = input$Othergenes3
        TCGAname = rep("TCGA", 17814)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:585) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }
      
      
      #
      # if (input$Endpoint=="Time-to-recurrence"){
      #   TCGA$Event=TCGA$recurrence
      #   TCGA$Time_To_Event=TCGA$time_to_recurrence
      # }
      # else if (input$Endpoint=="Overall survival"){
      #   TCGA$Event=TCGA$Event
      #   TCGA$Time_To_Event=TCGA$Time_To_Event_TO_Event}
      dd <<- datadist(TCGA)
      options(datadist = 'dd')
      #survival estimates stratified on Race
      f.ia = cph(
        Surv(Time_To_Event, Event) ~ rcs(Age, 3) + strat (Target_Gene_Stratified),
        x = TRUE ,
        y = TRUE ,
        surv = TRUE,
        data = TCGA
      )
      p1 = Predict(f.ia , Age , Target_Gene_Stratified , time =
                     365)
      p1 = ggplot(p1, aes(col = Target_Gene_Stratified)) +
        geom_line(aes(group = factor(Target_Gene_Stratified)), size = 2) +
        labs(title = paste(
          "Predicted Survival Probability for 1 Year",
          names(TCGA$target_gene)
        ))
      p2 = Predict(f.ia , Age , Target_Gene_Stratified , time =
                     365 * 2)
      p2 = ggplot(p2, aes(col = Target_Gene_Stratified)) +
        geom_line(aes(group = factor(Target_Gene_Stratified)), size = 2) +
        labs(title = paste(
          "Predicted Survival Probability for 2 Years",
          names(TCGA$target_gene)
        ))
      grid.arrange(p1, p2, nrow = 1, ncol = 3)
    })
  })
  #output$ui<-renderUI({
  observe({
    final_formu = sfreact()$final_formu
    
    a = sfreact()$a
    varlist = sfreact()$varlist
    #print (varlist)
    print (length(varlist[[1]]))
    # for (i in 1:length(varlist[[1]])){
    updateTextInput(session,
                    "inText",
                    label = paste(
                      "Please enter a value for:",
                      paste(unlist(varlist[[1]]), collapse = ","),
                      " separated
                      by space"
                    ))
    
    #value = as.character(varlist[[1]][i])
    #}
    
    #updateTextInput(session, "inText","Put your numbers")
    #print (length(varlist[[1]]))
    # for (i in 1:length(varlist[[1]])){
    # updateTextInput(session, paste0("inText", i), label = "Variable", value = as.character(varlist[[1]][i]))}
    #                                          #                   label = paste("Please enter a value:",varlist[[1]][i])))
    # Can also set the label, this time for input$inText2
    #   for (i in 1:length(varlist)){
    #   updateTextInput(session, paste("inText",i),
    #                   label = paste("Please enter a value:",varlist[[1]][i]))
    # }
    # df=setNames(data.frame(matrix(ncol = length(x), nrow = 0)), unlist(x))
    # df[1,]=unlist(y)})
    #})
  })
  
  
  output$Predplot2 <- renderPlot({
    #if input$Model==With Restricted Cubic Spline
    #print (TCGA$Time_To_Event)
    input$coxcalc
    isolate({
      final_formu = sfreact()$final_formu
      print (final_formu)
      # rcslist = sfreact()[[2]]
      #print (rcslist)
      a = sfreact()$a
      varlist = sfreact()$varlist
      #formu=paste("Surv(Time_To_Event, Event)~", unlist(varlist),collapse = "+",sep = "")
      #print (formu)
      # TCGA$target_gene = calcultatecox()
      TCGA = calcultatecox()[[1]]
      #print (TCGA[,1])
      TCGA$target_gene = calcultatecox()[[2]]
      
      if (input$Database == "TCGA GBM RNA-Seq") {
        targetname = input$Multigene2
        Othergenes = input$Othergenes2
        TCGAname = rep("TCGA", 20028)
        medianSelected = median(TCGA$target_gene)
        
        for (i in 1:172) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      } else if (input$Database == "TCGA GBM U133a") {
        targetname = input$Multigene
        Othergenes = input$Othergenes
        TCGAname = rep("TCGA", 12042)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:539) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }  else if (input$Database == "TCGA GBM Agilent") {
        targetname = input$Multigene3
        Othergenes = input$Othergenes3
        TCGAname = rep("TCGA", 17814)
        
        medianSelected = median(TCGA$target_gene)
        for (i in 1:585) {
          if (TCGA$target_gene[i] < medianSelected) {
            TCGA$Target_Gene_Stratified[i] = "LOW"
          }
          else {
            TCGA$Target_Gene_Stratified[i] = "HIGH"
          }
        }
      }
      
      datalist = strsplit(as.character(req(input$inText)), " ")
      #print (datalist)
      mat.or.vec(3, 2)
      newdataframe = data.frame(mat.or.vec(1, length(varlist[[1]])))
      
      # newdataframe= data.frame(matrix(,nrow=1,ncol=length(varlist[[1]])))
      #print (length(varlist))
      # print (length(varlist[[1]]))
      for (i in (1:length(varlist[[1]]))) {
        colnames(newdataframe)[i] <- varlist[[1]][i]
      }
      # for datalist, convert some character/string to numeric
      for (i in 1:length(datalist[[1]])) {
        if (varlist[[1]][i] %in% c(
          "KPS",
          "Gender",
          "GeneExp_Subtype",
          "G_Cimp_Status",
          "Initial_Pathologic",
          "Radiation_Therapy",
          "Race",
          "Target_Gene_Stratified"
        ))
        {
          newdataframe[1, i] = toupper(datalist[[1]][i])
        }
        else {
          newdataframe[1, i] = as.numeric(datalist[[1]][i])
        }
      }
      #print (length(datalist[[1]]))
      #print (str(newdataframe[1,5]))
      # for (i in 1:length(datalist[[1]])){
      # newdataframe[1,i]=unlist(datalist[[1]][i])
      # }
      #strsplit(as.character(input$inText2), " ")
      #print (typeof(newdataframe[1,5]))
      print (newdataframe)
      #print (varlist[[1]][1])
      #rcs or not:
      #if if (length(grep ("\\<rcs\\>", formula))>0){
      f = plot(
        survfit(cph(
          as.formula(final_formu), TCGA, x = TRUE, y = TRUE
        ), newdataframe),
        col = c("red", "green", "green"),
        xscale = 365.25,
        xlab = "Years",
        ylab = "Survival Probability"
      )
      legend(
        100,
        .2,
        c("Predicted Survival Probability", "95% Confidence Interval"),
        lty = c(1:2),
        col = c("red", "green")
      )
      
    })
  })
  #   input$inText2<-"Hello"
  # )
  ###################################### Module-- Univariate Analysis#########################################################
  ## 20180131: Two plots, so we need to have a reactive function
  calcultate_uni <-eventReactive(input$Uni_cox,
                                isolate({if (input$Databases == "TCGA GBM RNA-Seq") {
                                  TCGA = TCGA_GBM_rna_seq_overall
                                  target = TCGA[, input$Unigene2]
                                  name = input$Unigene2
                                } else if (input$Database == "TCGA GBM U133a") {
                                  TCGA = TCGA_GBM_u133a_overall
                                  target = TCGA[, input$Unigene]
                                  name = input$Unigene
                                
                                } else if (input$Database == "TCGA GBM Agilent") {
                                  TCGA = TCGA_GBM_agilent_overall
                                  target = TCGA[, input$Unigene3]
                                  name = input$Unigene3
                                }
                                  #print (TCGA[,1])
                                  return(list(TCGA, target,name))
                                }))
  
  recursive_tfit <- reactive({
    input$Uni_cox
    isolate({TCGA = calcultate_uni()[[1]]
    #print (TCGA[,1])
    TCGA$target = calcultate_uni()[[2]]
    tfit = rpart(formula = Surv(Time_To_Event, event = Event) ~ target, data = TCGA, control=rpart.control(minsplit=30, cp=0.01))
    return (tfit)
    })
  })
#   recursive_tfit <- reactive({
#     TCGA = correlation_database()
#     selected_gene=TCGA[,input$Unigene]
#     tfit = rpart(formula = Surv(Time_To_Event, event = Event) ~ selected_gene, data = TCGA, control=rpart.control(minsplit=30, cp=0.01))
#   })
  output$recurplot<-renderPlot({
    tfit2<-as.party(recursive_tfit())
    plot(tfit2,
         main="Conditional Inference Tree Plot")
  })
  output$uniplot <- renderPlot({
    TCGA = calcultate_uni()[[1]]
    #print (TCGA[,1])
    target_gene2 = calcultate_uni()[[2]]
    name = calcultate_uni()[[3]]
    tfit=recursive_tfit()
    #tfit = rpart(formula = Surv(TIME_TO_EVENT, event = EVENT) ~ selected_gene, data = TCGA, control=rpart.control(minsplit=30, cp=0.01))
    #tfit2 <- as.party(tfit)
    #plot(tfit2,
    #main="Conditional Inference Tree Plot")
    splitpoint=tfit$splits[1,4]
    TCGA[,"Split"]=NA
    TCGA$Split[which(target_gene2>splitpoint)]="High"
    TCGA$Split[which(target_gene2<=splitpoint)]="Low"
    #TCGA$Split=as.factor(TCGA$Split)
    print(TCGA$Split)
    km.by.split<- survfit(Surv(Time_To_Event,  Event) ~ Split, data = TCGA)
    #print(km.by.split)
    ggsurv=ggsurvplot(km.by.split, TCGA, # Change legends: title & labels
                      legend.title = name,
                      legend.labs = c("LOW", "HIGH"),
                      # Add p-value and tervals
                      pval = TRUE,
                      pval.method = TRUE,
                      log.rank.weights="n",
                      #pval.method.size=TRUE,
                      #conf.int = TRUE,
                      # Add risk table
                      risk.table = TRUE,
                      surv.scale="percent",
                      tables.height = 0.2,
                      tables.theme = theme_cleantable(),
                      # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                      # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                      palette = c( "#2E9FDF","#E7B800"),
                      font.x = 12, font.y = 12, font.main = 14, ylab = "% Surviving",
                      title = "Kaplan Meier Survival Estimates for selected gene high vs low")
    ggsurv$plot=ggsurv$plot+theme(plot.title = element_text(hjust = 0.5))
    print(ggsurv)
    
  })
})

