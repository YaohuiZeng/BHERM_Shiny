# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.


options(shiny.maxRequestSize = 9*1024^2)
source("utilities.R")
source("main.R")
source("post_analysis.R")
source("utilities_jags.R")
source("main_jags.R")
source("post_analysis_jags.R")

shinyServer(function(input, output, session) {

  dat <- reactive({
    infile1 <- input$file1
    if (is.null(infile1)) {
      return(NULL)
    }
    read.csv(infile1$datapath, header = T, sep = ",")
  })
  
  observe({
    data <- dat()
    choices <- as.list(c("NA", names(data)))

    updateSelectInput(session, "id", choices = choices)
    updateSelectInput(session, "trt", choices = choices)
    updateSelectInput(session, "dose", choices = choices)
    updateSelectInput(session, "safe", choices = choices)
    updateSelectInput(session, "safe.pk", choices = choices)
    updateSelectInput(session, "eff", choices = choices)
    updateSelectInput(session, "eff.pk", choices = choices)
    updateSelectInput(session, "eff.time", choices = choices)
  })

  
  ## model fitting
  res <- reactiveValues(main.fit = NULL, plots = NULL, 
                        n.plots = NULL, run.stan = TRUE)
  
  ## running RStan
  observeEvent(input$run, {

    res$main.fit <- NULL
    res$plots <- NULL
    res$n.plots <- 0
    res$run.stan <- TRUE
    
    data.file <- input$file1$name
    prior.file <- input$file2$name
    safe.model <- input$safe.model
    eff.model <- input$eff.model
    chains <- input$chains
    warmup <- input$warmup
    iter <- input$iter
    thin <- input$thin
    seed <- input$seed
  
    model.case <- get_model_case(safe.model = safe.model, eff.model = eff.model) 

    if (model.case == 1) { 
      # Case 1: two responses, no cox
      varnames <- c(input$id, input$trt, input$dose, input$safe,
                    input$safe.pk, input$eff, input$eff.pk)
      plot.funs <- list(plot_dose_pk_EFF = plot_dose_pk_EFF,
                        plot_predict_EFF = plot_predict_EFF,
                        plot_dose_pk_SAF = plot_dose_pk_SAF,
                        plot_predict_SAF = plot_predict_SAF)
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred),
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))

    } else if (model.case == 2) {
      # Case 2: only safety model, CANNOT be cox
      varnames <- c(input$id, input$trt, input$dose, input$safe,
                    input$safe.pk)
      plot.funs <- list(plot_dose_pk_SAF = plot_dose_pk_SAF,
                        plot_predict_SAF = plot_predict_SAF)
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))
      
    } else if (model.case == 3) {
      # Case 3: only efficacy model, not Cox.
      varnames <- c(input$id, input$trt, input$dose, input$eff, input$eff.pk)
      plot.funs <- list(plot_dose_pk_EFF = plot_dose_pk_EFF,
                        plot_predict_EFF = plot_predict_EFF)
      
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))
      
    } else if (model.case == 4) {
      # Case 4: only Efficacy cox model
      varnames <- c(input$id, input$trt, input$dose, 
                    input$eff, input$eff.time, input$eff.pk)
      plot.funs <- list(plot_dose_pk_EFF = plot_dose_pk_EFF,
                        plot_mean_survival = plot_mean_survival)
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))
#       plot.names <- c("dose.pk.eff", "pred.eff", "pred.median", "pred.surv")
    } 
    
    if (any(is.na(varnames))) {
      stop("Missing some variables.")
    }
    
    cat("\nModel fitting - Start: ", format(Sys.time()), '\n')

#     log <- capture.output(
#     if(input$run) {
      # run Rstan
    if (model.case != 4) {
      fit <- main.shiny.fit(file = data.file, varnames = varnames, 
                            safe.model = safe.model,
                            eff.model = eff.model, 
                            pred.doses = 1100,
                            file.prior = prior.file, chains = chains,
                            warmup = warmup, iter = iter, thin = thin,
                            seed = seed)
    } else {
      fit <- main.shiny.fit(file = data.file, varnames = varnames, 
                            safe.model = safe.model,
                            eff.model = eff.model, 
                            file.prior = prior.file, chains = chains,
                            warmup = warmup, iter = iter, thin = thin,
                            seed = seed)
    }

    cat("\nModel fitting - End: ", format(Sys.time()), '\n')
    
    res$main.fit <- fit
    res$n.plots <- length(plot.funs)

    if (!is.null(res$main.fit)) {
      res$plots <- lapply(1:length(plot.funs), function(i){
        plot.pars <- plot.label[[i]]
        plot.pars[['fit']] <- fit
        do.call(plot.funs[[i]], plot.pars)
      })
      
      lapply(1:res$n.plots, function(i){
        output[[paste("plot", i, sep="") ]] <- renderPlot({
          print(res$plots[[i]])
        })
      })
      
    }

  })


  observeEvent(input$run.jags, {
    
    res$main.fit <- NULL
    res$plots <- NULL
    res$n.plots <- 0
    res$run.stan <- FALSE
    
    data.file <- input$file1$name
    prior.file <- input$file2$name
    safe.model <- input$safe.model
    eff.model <- input$eff.model
    chains <- input$chains
    warmup <- input$warmup
    iter <- input$iter
    thin <- input$thin
    seed <- input$seed
    
    model.case <- get_model_case_JAGS(safe.model = safe.model, eff.model = eff.model) 
    
    PKSAFcase<-ifelse(safe.model=='Logistic','SAFLogistic',ifelse(safe.model=='Poisson','SAFPoisson',ifelse(safe.model=='NegBin','SAFNegBin',ifelse(safe.model=='Emax','SAFEmax','SAFNULL'))))
    PKEFFcase<-ifelse(eff.model=='Logistic','EFFLogistic',ifelse(eff.model=='Poisson','EFFPoisson',ifelse(eff.model=='NegBin','EFFNegBin',ifelse(eff.model=='Emax','EFFEmax','EFFNULL'))))
    PKSAFEFFcases<-c(PKSAFcase,PKEFFcase) 
    
    if (model.case == 1) { 
      # Case 1: two responses, no cox
      varnames <- c(input$id, input$trt, input$dose, input$safe,
                    input$safe.pk, input$eff, input$eff.pk)
      plot.funs <- list(plot_dose_pk_EFF_JAGS = plot_dose_pk_EFF_JAGS,
                        plot_predict_EFF_JAGS = plot_predict_EFF_JAGS,
                        plot_dose_pk_SAF_JAGS = plot_dose_pk_SAF_JAGS,
                        plot_predict_SAF_JAGS = plot_predict_SAF_JAGS)
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred),
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))
      
    } else if (model.case == 2) {
      # Case 2: only safety model, CANNOT be cox
      varnames <- c(input$id, input$trt, input$dose, input$safe,
                    input$safe.pk)
      plot.funs <- list(plot_dose_pk_SAF_JAGS = plot_dose_pk_SAF_JAGS,
                        plot_predict_SAF_JAGS = plot_predict_SAF_JAGS)
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))
      
    } else if (model.case == 3) {
      # Case 3: only efficacy model, not Cox.
      varnames <- c(input$id, input$trt, input$dose, input$eff, input$eff.pk)
      plot.funs <- list(plot_dose_pk_EFF_JAGS = plot_dose_pk_EFF_JAGS,
                        plot_predict_EFF_JAGS = plot_predict_EFF_JAGS)
      
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))
      
    } else if (model.case == 4) {
      # Case 4: only Efficacy cox model
      varnames <- c(input$id, input$trt, input$dose, 
                    input$eff, input$eff.time, input$eff.pk)
      plot.funs <- list(plot_dose_pk_EFF_JAGS = plot_dose_pk_EFF_JAGS,
                        plot_mean_survival_JAGS = plot_mean_survival_JAGS)
      plot.label <- list(
        list(xlab = input$xlab_pk, ylab = input$ylab_pk),
        list(xlab = input$xlab_pred, ylab = input$ylab_pred))
      #       plot.names <- c("dose.pk.eff", "pred.eff", "pred.median", "pred.surv")
    } 
    
    if (any(is.na(varnames))) {
      stop("Missing some variables.")
    }
    
    cat("\nModel fitting - Start: ", format(Sys.time()), '\n')
    
    #     log <- capture.output(
    #     if(input$run) {
    # run Rstan
    if (model.case != 4) {
      fit <- main_JAGS_shiny(file = data.file, varnames = varnames, 
                            safe.model = safe.model,
                            eff.model = eff.model, 
                            PKSAFEFFcases = PKSAFEFFcases,
#                             pred.doses = 1100,
                            file.prior = prior.file, chains = chains,
                            warmup = warmup, iter = iter, thin = thin,
                            seed = seed)
    } else {
      fit <- main_JAGS_shiny(file = data.file, varnames = varnames, 
                            safe.model = safe.model,
                            eff.model = eff.model, 
                            file.prior = prior.file, chains = chains,
                            warmup = warmup, iter = iter, thin = thin,
                            seed = seed)
    }
    
    cat("\nModel fitting - End: ", format(Sys.time()), '\n')
    
    res$main.fit <- fit
    res$n.plots <- length(plot.funs)
    
    if (!is.null(res$main.fit)) {
      res$plots <- lapply(1:length(plot.funs), function(i){
        plot.pars <- plot.label[[i]]
        plot.pars[['fit']] <- fit
        do.call(plot.funs[[i]], plot.pars)
      })
      
      lapply(1:res$n.plots, function(i){
        output[[paste("plot", i, sep="") ]] <- renderPlot({
          print(res$plots[[i]])
        })
      })
      
    }
    
  })

  #   values <- reactiveValues()
  #   log <- eventReactive(input$run, {
  #     fit <- res$main.fit
  #     
  #     get_summary_stat(fit$main.fit)
  # #     values[['log']] <- capture.output(res <- main.shiny())
  #   })


  ##### Create divs######
  output$pred.plots <- renderUI({
    if (!is.null(res$main.fit)) {
      plot_output_list <- lapply(1:res$n.plots, function(i) {
        plotname <- paste("plot", i, sep="")
        plotOutput(plotname, height = 400, width = 600)
      })   
      do.call(tagList, plot_output_list)
    } else {
      return(NULL)
    }
  })

  # observe({
  #   lapply(1:res$n.plots, function(i){
  #     output[[paste("plot", i, sep="") ]] <- renderPlot({
  #       print(res$plots[[i]])
  #     })
  #   })
  # })

  output$data <- renderTable({

    data <- dat()
    data
   
  })

  


  output$message <- renderPrint({
    cat("------------------------------------------------------------\n")
    cat("BHERM Shiny Instruction: \n\n")
    cat("\tStep 1: Choose models of Exposure-Response for EFFF and/or SAF.\n")
    cat("\tStep 2: Upload CSV file for raw data.\n")
    cat("\tStep 3: Select variables.\n")
    cat("\tStep 4: Upload R file for prior specifications.\n")
    cat("\tStep 5: (optional) Specify sampling parameters and others.\n\n")
    cat("Refer to Manual for detailed model specifications, parameters, and user guide.")
    cat("\n------------------------------------------------------------\n\n")
    if (input$safe.model == "NULL") {
      cat("SAF model selected: NA\n")
    } else {
      cat("SAF model selected: ", input$safe.model, "\n")
    }
    if (input$eff.model == 'NULL') {
      cat("EFF model selected: NA\n\n")
    } else {
      cat("EFF model selected: ", input$eff.model, "\n\n")
    }

  })

  output$console <- renderPrint({
    fit <- res$main.fit
    run.stan <- res$run.stan
    pars.model <- fit$pars.model
    
    names_summary_stat <- colnames(fit$jags.fit$sims.array[,1,]) 
    
    npars <- length(names_summary_stat)
    
    npars.model <- unlist(lapply(pars.model, function(x) {
      grep(pattern = x, names_summary_stat)
    }))
    
    
    cat("------------------------------------------------------------\n")
    cat("Posterior summary statistics: ")
    cat("\n------------------------------------------------------------\n\n")
#     if (!is.null(res$main.fit)) {
    if (run.stan) {
      print(fit$stan.fit, pars = fit$pars.model, digits_summary = 4)
    } else {
      ## YZ: need change to just ouput model parameters.
      print(fit$jags.fit, pars = fit$pars.model, digits_summary = 4)
    }
      
#     }
#     log()
#     return(print(values[['log']]))
  })

  output$trace.plot <- renderPlot({
    if (!is.null(res$main.fit)) {
      fit <- res$main.fit
      pars.model <- fit$pars.model
      pars.model <- c(pars.model, "lp__")
      
      ## YZ: need add option to output traceplots for JAGS
      stan_trace(fit$stan.fit, pars = pars.model) +
        ggtitle("Trace plot")
    }
  })

  output$dens.plot <- renderPlot({
    if (!is.null(res$main.fit)) {
      fit <- res$main.fit
      pars.model <- fit$pars.model
      pars.model <- c(pars.model, "lp__")
      
      ## YZ: need add option to output density plots for JAGS
      stan_dens(fit$stan.fit, pars = pars.model) +
        ggtitle("Density plot")
    }
  })



})

