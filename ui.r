shinyUI(  
  
  navbarPage(title = strong("BrainGeneExpress"), windowTitle = "", 
             fluid = TRUE, 
             
             source("tabs/CorrelationTab.r", local = TRUE)$value,
             source("tabs/MultivariableTab.r", local = TRUE)$value,
             source("tabs/UnivariableTab.r", local = TRUE)$value
             
  )
)
