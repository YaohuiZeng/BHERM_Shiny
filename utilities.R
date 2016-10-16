## all utility functions

## need install a package - ggrepel
# devtools::install_github("slowkow/ggrepel")

## load packages
#--------------------------------------------------------------------------
if(!require(ggrepel)) {
  if (!require(devtools)) {
    install.packages("devtools")
  }
  devtools::install_github("slowkow/ggrepel")
}
require(ggrepel)
require(gridExtra)
require(grid)
require(ggplot2)
require(rstan)
require(coda)

load("data/eff.cox.stan.fit0.RData")
load("data/eff.only.stan.fit0.RData")
load("data/safe.eff.stan.fit0.RData")
load("data/safe.only.stan.fit0.RData")
load("data/obs.prob.RData")
# require.packages <- c("ggplot2", "gridExtra", "grid", "rstan", "coda")
# lapply(require.packages, function(x) {
#   if (!require(x)) {
#     install.packages(x)
#     require(x)
#   }
# })
#--------------------------------------------------------------------------

#==========================================================================
## utilities functions
#==========================================================================
## function to convert stan to coda
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[, x, ])))
}

# data in source file to list
source2list <- function(file){ 
  source(file, local = TRUE) 
  d <- mget(ls()) 
  d$file <- NULL 
  d 
} 


## interpolate survival time corresponding to quart 
interp_quart <- function(x, time_seq, quart) {
  id_upper <- which(x > quart)
  id_lower <- which(x < quart)
  id_quart <- which(x == quart)
  
  if (length(id_quart)) {
    return(time_seq[id_quart])
  } else if (!length(id_upper) || !length(id_lower)) {
    return(NA)
  } else {
    id_up <- id_upper[length(id_upper)]
    id_low <- id_lower[1]
    t <- (time_seq[id_up] - time_seq[id_low]) / (x[id_up] - x[id_low])  *
      (quart - x[id_low]) + time_seq[id_low]
    
    return (t)
  }
}


summary_data_cox <- function(dat, rescale) {
  # input: dat in function summary_data
  # output: summary data for cox model

  # EFF is censor indicator 
  names(dat) <- c('ID', 'TRT', 'DOSE', 'EFF', 'EFF_TIME', 'PK_EFF')
  na.id <- which(is.na(dat$EFF_TIME))
  dat <- dat[-na.id, ] # remove obs with missing follow-up time
  
  N <- nrow(dat)
  TRT <- levels(dat[, 2])
  NdA <- length(unique(TRT))
  dA <- unique(sort(dat[, 3]))
  dAn <- min(dA) # reference dose
  
  # if missing censor, set to be 1
  dat$EFF <- ifelse(is.na(dat$EFF), 1, dat$EFF)
  # call function to summarize efficacy data
  dat.EFF <- summary_data_single(dat = dat[, -5], rescale = rescale, is.SAF = FALSE)
  # drop EFFprob_dose_cat, NEFF, EFF, not need for Cox model
  dat.EFF$EFFprob_dose_cat <- NULL
  dat.EFF$EFF <- NULL
  dat.EFF$NEFF <- NULL
  
  event.t <- with(dat, sort(unique(EFF_TIME[EFF == 0])))
  obs_t <- dat$EFF_TIME
  t <- c(event.t, max(dat$EFF_TIME, na.rm = TRUE))
  NT <- length(event.t)

  eps <- 1.0e-10
  fail <- 1- dat$EFF
  
  # Define Y and dN
  Y <- matrix(NA, NdA, NT)
  dN <- matrix(NA, NdA, NT)
  Ytemp <- matrix(NA, N, NT)
  dNtemp <- matrix(NA, N, NT)
  
  # counting process per patient
  for (p in 1:N){ 
    for (q in 1:NT){             
      Ytemp[p, q] <- ifelse(obs_t[p] - t[q] + eps < 0, 0, 1)
      dNtemp[p, q] <- Ytemp[p, q] * ifelse(t[q + 1] - obs_t[p] - eps < 0, 0, 1) * fail[p]
    }   
  }
  # counting process per dose
  for (i in 1:NdA){
    for (j in 1:NT){
      sub<-NULL
      sub<-which(dat$DOSE == dA[i])
      Y[i,j]<-sum(Ytemp[sub, j])
      dN[i,j]<-sum(dNtemp[sub, j])
    }
  }
  
  data <- c(list(TRT = TRT, NT = NT, NdA = NdA, dA = dA, dAn = dAn, t = t,
                 Y = Y, dN = dN), dat.EFF)
  data
  
}


## prepare stan data for Efficacy Cox model
prepare_stan_data_cox <- function(summary.data, 
                                  file.prior = NULL,
                                  pred.doses = NULL, 
                                  pred.time.max = NULL, # max prediction time
                                  increment.dose = 10,
                                  increment.time = 10,
                                  pred.NPK = 3) {
  # get prior data from file
  if (is.null(file.prior)) {
    file.prior <- 'default_priors.R'
  }
  
  dat.prior.input <- source2list(file.prior)
  dat.prior <- list()
  
  # prior for DOSE-PK model
  dat.prior$mulogphiA_EFF <- dat.prior.input$mulogphiA_EFF
  dat.prior$Sigma_logphiA_EFF <- dat.prior.input$Sigma_logphiA_EFF
  dat.prior$mulogsigmaA_EFF <- dat.prior.input$mulogsigmaA_EFF
  dat.prior$sigma_logsigmaA_EFF <- dat.prior.input$sigma_logsigmaA_EFF
  
  # prior for Cox model
  dat.prior$muEFF_cox <- dat.prior.input$muEFF_cox
  dat.prior$SigmaEFF_cox <- dat.prior.input$SigmaEFF_cox
  dat.prior$c <- dat.prior.input$c
  dat.prior$r <- dat.prior.input$r
  
  # posterior sampling over a finer grid of doses within range of DosesA
  if (is.null(pred.doses)) {
    DosesA <- summary.data$dA
  } else {
    DosesA <- pred.doses
  }
  if (is.null(pred.time.max)) {
    pred.time.max <- max(summary.data$t) * 2
  } else {
    pred.time.max <- pred.time.max
  }
  ## TODO: could be 0 if placebo!
  DoseRefA <- min(DosesA)
  # posterior sampling over a finer grid of doses within range of DosesA
  seq_dosesA <- seq(from = min(DosesA), to = max(DosesA), by = increment.dose)
  Dose_len <- length(seq_dosesA)
#   indice_DosesA <- which(seq_dosesA %in% DosesA)
  seq_NPK <- rep(pred.NPK, Dose_len)
  
  # posterior sampling over a finer grid of time within range of c(0, pred.time.max)
  seq_t_pred <- seq(0, pred.time.max, by = increment.time)
  if (max(seq_t_pred) < pred.time.max) {
    seq_t_pred <- c(seq_t_pred, pred.time.max)
  }
  T_len <- length(seq_t_pred) - 1
  
  summary.dat <- summary.data[!(names(summary.data) %in% c("TRT", "rescale.pk.safe", "rescale.pk.eff"))]
  
  stan.data <- c(summary.dat, dat.prior, 
                 list(Dose_len = Dose_len, seq_dosesA = seq_dosesA, 
#                       indice_DosesA = indice_DosesA, DosesA = DosesA, 
                      seq_NPK = seq_NPK, 
                      DoseRefA = DoseRefA, 
                      seq_t_pred = seq_t_pred, T_len = T_len))
  
  stan.data
}


summary_data_single <- function(dat, rescale, is.SAF) {
  if (is.SAF) {
#     names(dat) <- c('ID', 'TRT', 'DOSE', 'SAF', 'PK_SAF')
    dat$logPK_SAF <- log(dat$PK_SAF)
    
    if(is.null(rescale)) {
      rescale.pk.safe <- exp(mean(dat$logPK_SAF, na.rm = TRUE))
      dat$logPK_SAF <- dat$logPK_SAF - mean(dat$logPK_SAF, na.rm = TRUE)
    } else {
      rescale.pk.safe <- rescale
      dat$logPK_SAF <- dat$logPK_SAF - log(rescale)
    }

    # cut-off for SAF regions
    SAFprob_dose_cat <- c(0, 0.16, 0.33, 1)
    
    ## aggregate function omits NA by default
    NSAF <- aggregate(SAF ~ TRT, data = dat, length)$SAF
    SAF <- aggregate(SAF ~ TRT, data = dat, function(x) sum(x == 1))$SAF
    NPK_SAF <- aggregate(PK_SAF ~ TRT, data = dat, length)$PK_SAF
    
    logX_SAF <- aggregate(logPK_SAF ~ TRT, data = dat, 
                          mean, na.rm = T)$logPK_SAF 
    samplevar_SAF <- aggregate(logPK_SAF ~ TRT, data = dat, 
                               var, na.rm = T)$logPK_SAF  
    Y_SAF <- samplevar_SAF*(NPK_SAF-1)
    Y_SAF <- ifelse(is.na(Y_SAF), 1e10, Y_SAF) # remove NA
    
    data <- list(
      NSAF = NSAF, SAF = SAF, NPK_SAF = NPK_SAF, 
      logX_SAF = logX_SAF, Y_SAF = Y_SAF,
      rescale.pk.safe = rescale.pk.safe,
      SAFprob_dose_cat = SAFprob_dose_cat
    )
    return(data)
  } else {
#     names(dat) <- c('ID', 'TRT', 'DOSE', 'EFF', 'PK_EFF')
    dat$logPK_EFF <- log(dat$PK_EFF)
    if(is.null(rescale)) {
      rescale.pk.eff <- exp(mean(dat$logPK_EFF, na.rm = TRUE))
      dat$logPK_EFF <- dat$logPK_EFF - mean(dat$logPK_EFF, na.rm = TRUE)
    } else {
      rescale.pk.eff <- rescale
      dat$logPK_EFF <- dat$logPK_EFF - log(rescale)
    }
    
    # cut-off for EFF regions
    EFFprob_dose_cat <- c(0, 0.4, 1)
    
    ## aggregate function omits NA by default
    NEFF <- aggregate(EFF ~ TRT, data = dat,length)$EFF
    EFF <- aggregate(EFF ~ TRT, data = dat, function(x) sum(x == 1))$EFF
    NPK_EFF <- aggregate(PK_EFF ~ TRT, data = dat, length)$PK_EFF
    
    logX_EFF <- aggregate(logPK_EFF ~ TRT, data = dat, 
                          mean, na.rm = T)$logPK_EFF 
    samplevar_EFF <- aggregate(logPK_EFF ~ TRT, data = dat, 
                               var, na.rm = T)$logPK_EFF  
    Y_EFF <- samplevar_EFF*(NPK_EFF-1)
    Y_EFF <- ifelse(is.na(Y_EFF), 1e10, Y_EFF)
    
    data <- list( 
      NEFF = NEFF, EFF = EFF, NPK_EFF = NPK_EFF, 
      logX_EFF = logX_EFF, Y_EFF = Y_EFF,
      rescale.pk.eff = rescale.pk.eff,
      EFFprob_dose_cat = EFFprob_dose_cat
    )
    return(data)
  }
}


## get model specification scenario
get_model_case <- function(safe.model, eff.model) {
  if (safe.model == 'NULL' && eff.model == 'NULL') {
    stop("Must specify at least one model.")
  }
  
  if (safe.model %in% c("Logistic", "Poisson", "NegBin", "Emax") && 
        eff.model %in% c("Logistic", "Poisson", "NegBin", "Emax")) {
    # Case 1: two responses, no cox
    model.case <- 1
  } else if (safe.model %in% c("Logistic", "Poisson", "NegBin", "Emax") && eff.model == "NULL") {
    # Case 2: only safety model, CANNOT be cox
    model.case <- 2
  } else if (safe.model == 'NULL' && eff.model %in% c("Logistic", "Poisson", "NegBin", "Emax")) {
    # Case 3: only efficacy model, not Cox.
    model.case <- 3
  } else if (safe.model == 'NULL' && eff.model == "Cox") {
    # Case 4: only Efficacy cox model
    model.case <- 4
    
  } else if (safe.model %in% c("Logistic", "Poisson", "NegBin", "Emax") && eff.model == 'Cox') {
    # Case 5: Efficacy Cox model + Safety model
    model.case <- 5
    stop("The Safety model and Efficacy Cox model need fit separately. ")
  } else {
    stop("Invalid model option.")
  }
  
  model.case
}


## get summary data from file
summary_data <- function(file, varnames, 
                         safe.model = c("Logistic", "Poisson", "NegBin", "Emax", "NULL"), 
                         eff.model = c("Logistic", "Poisson", "NegBin", "Emax", "Cox", "NULL"),
                         rescale.pk.safe = NULL,
                         rescale.pk.eff = NULL) {
  
  safe.model <- match.arg(safe.model)
  eff.model <- match.arg(eff.model)
  
  model.case <- get_model_case(safe.model, eff.model)
  
  raw <- read.csv(file = file, header = TRUE, sep = ',')
  if (!all(varnames %in% names(raw))) {
    stop("varnames list doesn't match variable names in file!")
  }
  nvars <- length(varnames)
  dat <- raw[, varnames]
  
  TRT <- levels(dat[, 2])
  NdA <- length(unique(TRT))
  dA <- unique(sort(dat[, 3]))
  
  dAn <- min(dA) # reference dose
  dat.DOSE <- list(TRT = TRT, NdA = NdA, dA = dA, dAn = dAn)
  
  if (model.case == 1) {
    # Case 1: two responses, no cox
    if (nvars != 7) {
      stop("There must be 7 variables in 'varnames' for both Safety and Efficacy models. ")
    }

    names(dat) <- c('ID', 'TRT', 'DOSE', 'SAF', 'PK_SAF','EFF', 'PK_EFF')
    dat.SAF <- summary_data_single(dat = dat[, 1:5], rescale = rescale.pk.safe, is.SAF = TRUE)
    dat.EFF <- summary_data_single(dat = dat[, c(1:3, 6, 7)], 
                                   rescale = rescale.pk.eff, is.SAF = FALSE)
    data <- c(dat.DOSE, dat.SAF, dat.EFF)
    
    return(data)
    
  } else if (model.case == 2) {
    # Case 2: only safety model, CANNOT be cox
    if (nvars != 5) {
      stop("There must be 5 variables in 'varnames' for Safety-only model if not Cox. ")
    }

    names(dat) <- c('ID', 'TRT', 'DOSE', 'SAF', 'PK_SAF')
    dat.SAF <- summary_data_single(dat = dat, rescale = rescale.pk.safe, is.SAF = TRUE)
    data <- c(dat.DOSE, dat.SAF)
    
    return(data)
    
  } else if (model.case == 3) {
    # Case 3: only efficacy model, not Cox.
    if (nvars != 5) {
      stop("There must be 5 variables in 'varnames' for Efficacy model. ")
    }
    names(dat) <- c('ID', 'TRT', 'DOSE', 'EFF', 'PK_EFF')
    dat.EFF <- summary_data_single(dat = dat, rescale = rescale.pk.eff, is.SAF = FALSE)
    data <- c(dat.DOSE, dat.EFF)
    return(data)
    
  } else if (model.case == 4) {
    # Case 4: only Efficacy cox model
#     stop("Will do")

    if (nvars != 6) {
      stop("There must be 6 variables in 'varnames' for Efficacy Cox model.")
    }
    names(dat) <- c('ID', 'TRT', 'DOSE', 'EFF_CENSOR', 'EFF_TIME', 'PK_EFF')
    dat.EFF.cox <- summary_data_cox(dat = dat, rescale = rescale.pk.eff)
    return(dat.EFF.cox)
    ## TODO:
    ## names(dat) <- ...
    ## data <- get_summary_data(dat, rescale = rescale)
    ## return(data)

  } else if (model.case == 5) {
    # Case 5: Efficacy Cox model + Safety model
    stop("will do")
    
  } else {
    stop("Invalid model option.")
  }
}

## prepare for stan data
prepare_stan_data <- function(summary.data, 
                              safe.model, 
                              eff.model, 
                              file.prior = NULL,
                              pred.doses = NULL, 
                              increment = 10, 
                              pred.NPK = 3) {
  # get prior data from file
  if (is.null(file.prior)) {
    file.prior <- 'default_priors.R'
  }
  
  dat.prior.input <- source2list(file.prior)
  dat.prior <- list()
  
  if (safe.model == 'Logistic') {
    SAF_model <- 1
  } else if (safe.model == 'Poisson') {
    SAF_model <- 2
  } else if (safe.model == 'NegBin') {
    SAF_model <- 3
  } else if (safe.model == 'Emax') {
    SAF_model <- 4
  } else if (safe.model == 'Cox') { # surival
    SAF_model <- 5
    stop('Will do!')
  } else if (safe.model == 'NULL') {
    SAF_model <- NA
  } else {
    stop('Invalid Safety model option!')
  }
  
  if (safe.model %in% c('Logistic', 'Poisson', 'NegBin', 'Emax')) {
    dat.prior$mulogphiA_SAF <- dat.prior.input$mulogphiA_SAF
    dat.prior$Sigma_logphiA_SAF <- dat.prior.input$Sigma_logphiA_SAF
    dat.prior$mulogsigmaA_SAF <- dat.prior.input$mulogsigmaA_SAF
    dat.prior$sigma_logsigmaA_SAF <- dat.prior.input$sigma_logsigmaA_SAF
    dat.prior$alpha_logalpha1 <- dat.prior.input$alpha_logalpha1
    dat.prior$beta_logalpha1 <- dat.prior.input$beta_logalpha1
    dat.prior$SAF_model <- SAF_model
    
    if (safe.model == 'Emax') {
      dat.prior$muSAF <- dat.prior.input$muSAF_emax # 1-by-4 vector for Emax
      dat.prior$SigmaSAF <- dat.prior.input$SigmaSAF_emax # 4-by-4 matrix for Emax
    } else {
      dat.prior$muSAF <- dat.prior.input$muSAF # 1-by-2 vector for others
      dat.prior$SigmaSAF <- dat.prior.input$SigmaSAF # 2-by-2 matrix for others
    }
    dat.prior$npara_SAF <- length(dat.prior$muSAF)
  } else if (safe.model == 'Cox') {
    stop("will do")
  }
  
  if (eff.model == 'Logistic') {
    EFF_model <- 1
  } else if (eff.model == 'Poisson') {
    EFF_model <- 2
  } else if (eff.model == 'NegBin') {
    EFF_model <- 3
  } else if (eff.model == 'Emax') {
    EFF_model <- 4
  } else if (eff.model == 'Cox') { # surival
    EFF_model <- 5
    #     stop('Will do!')
  } else if (eff.model == 'NULL') {
    EFF_model <- NA
  } else {
    stop('Invalid Safety model option!')
  }
  
  if (eff.model %in% c('Logistic', 'Poisson', 'NegBin', 'Emax')) {
    dat.prior$mulogphiA_EFF <- dat.prior.input$mulogphiA_EFF
    dat.prior$Sigma_logphiA_EFF <- dat.prior.input$Sigma_logphiA_EFF
    dat.prior$mulogsigmaA_EFF <- dat.prior.input$mulogsigmaA_EFF
    dat.prior$sigma_logsigmaA_EFF <- dat.prior.input$sigma_logsigmaA_EFF
    dat.prior$alpha_logalpha2 <- dat.prior.input$alpha_logalpha2
    dat.prior$beta_logalpha2 <- dat.prior.input$beta_logalpha2
    dat.prior$EFF_model <- EFF_model
    if (eff.model == 'Emax') {
      dat.prior$muEFF <- dat.prior.input$muEFF_emax # 1-by-4 vector for Emax
      dat.prior$SigmaEFF <- dat.prior.input$SigmaEFF_emax # 4-by-4 matrix for Emax
    } else {
      dat.prior$muEFF <- dat.prior.input$muEFF
      dat.prior$SigmaEFF <- dat.prior.input$SigmaEFF
    }
    dat.prior$npara_EFF <- length(dat.prior$muEFF)
    
  } else if (eff.model == 'Cox') {
    #     stop("will do")
  }
  
  # posterior sampling over a finer grid of doses within range of DosesA
  DosesA <- summary.data$dA
  if(!is.null(pred.doses)) {
    DosesA[1] <- min(summary.data$dA)
    DosesA[length(summary.data$dA)] <- max(pred.doses)
  }
  
  DoseRefA <- min(DosesA)
  #   DoseRefA <- min(summary.data$dA)
  
  # posterior sampling over a finer grid of doses within range of DosesA
  seq_dosesA <- seq(from = min(DosesA), to = max(DosesA), by = increment)
  seq_len <- length(seq_dosesA)
  indice_DosesA <- which(seq_dosesA %in% DosesA)
  seq_NPK <- rep(pred.NPK, seq_len)
  
  summary.dat <- summary.data[!(names(summary.data) %in% c("TRT", "rescale.pk.safe", "rescale.pk.eff"))]
  
  stan.data <- c(summary.dat, dat.prior, 
                 list(seq_len = seq_len, seq_dosesA = seq_dosesA, 
                      indice_DosesA = indice_DosesA, 
                      seq_NPK = seq_NPK, 
                      DosesA = DosesA, 
                      DoseRefA = DoseRefA))

  
  stan.data
}

# -------------------------------------------
