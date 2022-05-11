shinyUI(  
  
  navbarPage(title = strong("BrainGeneExpress"), windowTitle = "", 
             fluid = TRUE, 
             
             source("tabs/CorrelationTab.r", local = TRUE)$value,
             source("tabs/CoxRegressionTab.r", local = TRUE)$value,
             source("tabs/SurvivalTreeTab.r", local = TRUE)$value
             
  )
)
