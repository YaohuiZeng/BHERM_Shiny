shinyUI(fluidPage(
  titlePanel("Bayesian Hierarchical Exposure Response Modeling (BHERM)"),

  sidebarLayout(
    sidebarPanel(
      # Copy the line below to make a select box 
      tags$h4("Choose models:"),

      fluidRow(
        column(6, selectInput("safe.model", label = h5("Safety model:"), 
                              choices = list("Logistic" = "Logistic", 
                                             "Poisson" = "Poisson", 
                                             "NegBin" = "NegBin",
                                             "Emax" = "Emax",
                                             "Cox" = "CoX",
                                             "NA" = "NULL"), 
                              selected = "NULL"
#                               ,width = "50%"
        )),
        column(6, selectInput("eff.model", label = h5("Efficacy model:"), 
                              choices = list("Logistic" = "Logistic", 
                                             "Poisson" = "Poisson", 
                                             "NegBin" = "NegBin",
                                             "Emax" = "Emax",
                                             "Cox" = "Cox",
                                             "NA" = "NULL"), 
                              selected = "NULL")
               )

        ),
      fileInput('file1', 'Load Data file (.csv):',
                accept = c('.csv')
      ),
      
#       wellPanel(
        tags$h4("Select variables:"),
        selectInput("id", label = h6("Subject ID:"), 'NA'),
        fluidRow(
          column(6, selectInput("trt", label = h6("Treatment:"), 'NA')),
          column(6, selectInput("dose", label = h6("Dose level:"), "NA"))
          ),

        fluidRow(
          column(6, selectInput("safe", label = h6("Safety endpoint:"), "NA")),
          column(6, selectInput("eff", label = h6("Efficacy endpoint:"), "NA"))
          ),
        fluidRow(
          column(6, selectInput("safe.pk", label = h6("Safety PK parameter:"), "NA")),
          column(6, selectInput("eff.pk", label = h6("Efficacy PK parameter:"), "NA"))
          ),
        fluidRow(
          column(6, selectInput("safe.time", label = h6("Safety follow-up time:"), "NA")),
          column(6, selectInput("eff.time", label = h6("Efficacy follow-up time:"), "NA"))
          ),
        
      fileInput('file2', 'Piror file (.R):',
                accept = c('.R')
      ),
#       wellPanel(
        tags$h4("Sampling parameters:"),
        fluidRow(
          column(6, numericInput("chains", label = h6("Chains:"), 4)),
          column(6, numericInput("warmup", label = h6("Warmup:"), 500))
          ),
        fluidRow(
          column(6, numericInput("iter", label = h6("Iterations:"), 1500)),
          column(6, numericInput("thin", label = h6("Thin:"), 1))
        ),

        numericInput("seed", label = h6("Seed:"), 2016
#                      as.numeric(gsub("-", "", Sys.Date()))
                     ),

      tags$h4("Plot parameters:"),
      fluidRow(
        column(6, textInput("xlab_pk", label = h6("x label (Dose-PK):"), "Dose")),
        column(6, textInput("ylab_pk", label = h6("y lable (Dose-PK):"), "PK"))
      ),
      fluidRow(
        column(6, textInput("xlab_pred", label = h6("x label (PK-Response):"), "PK")),
        column(6, textInput("ylab_pred", label = h6("y lable (PK-Response):"), "Probability"))
      ),


    width = 4    
    ),

    mainPanel(

      fluidRow(
        column(3, 
               actionButton("run", "Run RStan", icon("paper-plane"),
                            style="color: #fff; background-color: #337ab7;
                     border-color: #2e6da4")),
        column(3, 
               actionButton("run.jags", "Run JAGS", icon("paper-plane"),
                            style="color: #fff; background-color: #337ab7;
                     border-color: #2e6da4")),
        column(3, 
               actionButton("pdf", "Manual", onclick = "window.open('BHERM.pdf')",
                            icon("hand-o-right"),
                            style="color: #fff; background-color: #337ab7;
                     border-color: #2e6da4"))
      ),
      tabsetPanel(type = 'tabs',
                  tabPanel("Console", 
                           fluidRow(column(12, verbatimTextOutput("message"))),
                           fluidRow(column(12, verbatimTextOutput("console")))),
                  tabPanel("Data", tableOutput('data')), 
                  tabPanel("Plot", uiOutput("pred.plots")),
                  tabPanel("Trace plot", 
                           plotOutput("trace.plot"),
                           plotOutput("dens.plot"))
                  )
    )
  )
))
