#' @title Performs the plasmode simulation
#' @description Creates 'plasmode' simulated datasets based on a given dataset when the outcome variable is time to event and exposure variable are binary. Plasmode simulation samples subjects with replacement from the observed data, uses subjects’ covariate data as is, and simulates exposure, outcome, or both.
#' @author Jessica M. Franklin, Younathan Abdia, and Shirley Wang
#' @param formulaOut An outcome model formula for estimating the hazard of outcome event.The functional form of the outcome model should be, Surv(data$time, data$event)~ Exposure + Confounders, where data is the dataset on which simulations are based, time is the follow-up time for the right-censored data and event is the status indicator. Exposure main effect must be first independent variable.
#' @param formulaCen An outcome model formula for estimating the hazard of censoring.The functional form of the outcome model should be, Surv(data$time, !data$event)~ Exposure + Confounders, where data is the dataset on which simulations are based, time is the follow-up time for the right-censored data and event is the status indicator. Exposure main effect must be first independent variable.
#' @param objectOut A fitted model object for the hazard of outcome.The functional form of the fitted model object should be of form coxph(Surv(data$time, data$event)~ Exposure + Confounders, data,x=TRUE), where coxph fits the Cox proportional hazard model, data is the dataset on which simulations are based, time is the follow-up time for the right-censored data and event is the status indicator. Exposure main effect must be first independent variable.
#' @param objectCen A fitted model object for the hazard of censoring.The functional form of the fitted model object should be of form coxph(Surv(data$time, !data$event)~ Exposure + Confounders, data,x=TRUE), where coxph fits the Cox proportional hazard model, data is the dataset on which simulations are based, time is the follow-up time for the right-censored data and event is the status indicator. Exposure main effect must be first independent variable.
#' @param formulaExp An exposure model formula containing the binary exposure on the left-hand side and potential confounders on the right-hand side. The functional form of the exposure model is, Exposure ~ Confounders.
#' @param objectExp A fitted model object for the exposure model.
#' @param data The dataset on which simulations are based.The data is required only when formulaOut, formulaCen or formulaExp or both are supplied to the argument.
#' @param idVar Name of the ID variable
#' @param effectOR The desired treatment effect odds ratio. By default effectOR = 1.
#' @param MMOut Multiplier of confounder effects on outcome on the log-scale. By default MMOut = 1 but one can specify a vector of length equivalent to the number of variables on the right-hand side of the outcome model.
#' @param MMExp Multiplier of confounder effects on exposure. By default MMExp = 1 but one can specify a vector of length equivalent to the number of variables on the right-hand side of the exposure model.
#' @param nsim Number of desired simulated datasets.
#' @param size Desired size of simulated datasets (i.e., # of individuals).
#' @param eventRate Desired average event rate. Default is the event rate in the observed data.
#' @param exposedPrev Desired average exposure rate. Default is the exposure prevalence in the observed data.
#' @details At least one of formulaOut, formulaCen, formulaExp, objectOut,objectCen, and objectExp must be specified, and which of these are specified will determine what gets simulated and how. If objectOut and objectCen or objectExp are specified, these objects are used as the base model for outcome and exposure simulation. If formulaOut and formulaCen or formulaExp are specified, then data should be given and base models are fit in the data using coxph with the given formulas. If formulaOut and formulaCen or objectOut and objectCen is specified, outcome will be simulated based on subjects’ observed exposure. If formulaExp or objectExp is specified, exposure will be simulated. And if models are specified for both outcome and exposure, both variables will be simulated with simulated outcome dependent on the simulated exposure.
#' @export
#' @import mgcv
#' @import glm2
#' @import nlme
#' @import survival
#' @importFrom stats as.formula coef fitted glm.control model.matrix plogis rbinom residuals rnorm runif uniroot var
#' @importFrom utils head tail
#' @return PlasmodeSur returns true beta coefficients used to generate the outcome and the exposure. PlasmodeSur also returns the data frame with the simulated data, including sampled IDs for each of nsim datasets along with simulated outcomes, exposure, or both.
#' \item{TrueOutBeta}{True beta coefficients used to generate the outcome.}
#' \item{TrueExpBeta}{True beta coefficients used to generate the exposure.}
#' \item{Sim_Data}{Plasmode simulated data, including sampled IDs for each of nsim datasets along with simulated outcomes, exposure, or both.}
#' @examples {
#' library(survival)
#' library(splines)
#' library(glm2)
#' ## Creating data set for simulation
#' lung <- lung[complete.cases(lung),]
#' lung$id <- 1:nrow(lung)
#' lung$meal.cal <- ifelse(lung$meal.cal > 1000, 1, 0)
#' lung$status <- lung$status - 1
#'
#' ## Formulas for estimating the hazard of outcome event, the hazard of censoring and exposure.
#'
#' form1<-Surv(lung$time, lung$status)~meal.cal+age+sex+ph.ecog+ph.karno
#' form2<-Surv(lung$time, !lung$status)~meal.cal+age+sex+ph.ecog+ph.karno
#' form3<- meal.cal~age+sex+ph.ecog+ph.karno
#'
#' Sur_Form1<-PlasmodeSur(formulaOut=form1,formulaCen=form2, objectOut=NULL, objectCen = NULL,
#'             formulaExp=form3,objectExp=NULL,data=lung,idVar="id",effectOR =1, MMOut=c(0.5,2,2,1,3),
#'             MMExp=c(2,2,2,2), nsim=3, size=nrow(lung), eventRate=NULL, exposedPrev=NULL)
#'
#' Sur_Form2<-PlasmodeSur(formulaOut=form1,formulaCen=form2, objectOut=NULL, objectCen = NULL,
#'             formulaExp=NULL,objectExp=NULL,data=lung,idVar="id",effectOR =1, MMOut=c(1,2,2,1,3),
#'             MMExp=c(1,1,1,1),nsim=3, size=nrow(lung), eventRate=NULL, exposedPrev=NULL)
#'
#' Sur_Form3<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=NULL, objectCen = NULL,
#'             formulaExp=form3,objectExp=NULL,data=lung,idVar="id",effectOR =1, MMOut=c(1,2,2,1,3),
#'             MMExp=c(1,1,1,1),nsim=3, size=nrow(lung), eventRate=NULL, exposedPrev=NULL)
#'
#' ## Objects for the hazard of the outcome event, hazard for censoring and the exposure.
#'
#' smod1 <- coxph(Surv(lung$time, lung$status)~meal.cal+age+sex+ph.ecog+ph.karno, data = lung,x=TRUE)
#' smod2 <- coxph(Surv(lung$time, !lung$status)~meal.cal+age+sex+ph.ecog+ph.karno, data = lung,x=TRUE)
#' pmod1<-glm2(meal.cal~age+sex+ph.ecog+ph.karno, data = lung,family = "binomial",
#'             control=glm.control(trace=TRUE))
#'
#' Sur_Obj1<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=smod1,objectCen = smod2,
#'             formulaExp=NULL,objectExp=pmod1,idVar=lung$id, effectOR =1, MMOut=c(1,2,2,1,3),
#'             MMExp=1, nsim=3,size=nrow(lung), eventRate=0.5, exposedPrev=NULL)
#'
#' Sur_Obj2<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=smod1,objectCen = smod2,
#'             formulaExp=NULL,objectExp=NULL,idVar=lung$id, effectOR =1.5, MMOut=c(1,2,2,1,3),
#'             MMExp=1, nsim=3,size=nrow(lung), eventRate=0.5, exposedPrev=NULL)
#'
#' Sur_Obj3<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=NULL,objectCen = NULL,
#'             formulaExp=NULL,objectExp=pmod1,idVar=lung$id,effectOR =1, MMOut=c(1,2,2,1,3),
#'             MMExp=1, nsim=3,size=nrow(lung), eventRate=0.5, exposedPrev=NULL)
#'}

PlasmodeSur<- function(formulaOut=NULL,formulaCen=NULL, objectOut=NULL,objectCen=NULL,formulaExp=NULL,objectExp=NULL,data,idVar,
                       effectOR =1, MMOut=1,MMExp=1, nsim, size, eventRate=NULL, exposedPrev=NULL)
{

  ## Code for simulating data when the outcome is Survival and data set is provided to estimated the outcome and exposure.

  if(is.null(formulaOut)==FALSE & is.null(formulaCen)==FALSE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE & is.null(objectCen)==TRUE & is.null(objectExp)==TRUE){
    n <- nrow(data)
    exposure<-all.vars(formulaOut)[4]
    sidx <- sapply(c(idVar),function(v) which(names(data) == v))
    names(data)[sidx] <- c("ID")

    x<-data

    ##estimate survival and censoring models
    smod <- coxph(formulaOut, x = TRUE, data = x)
    fit <- survfit(smod)
    s0 <- fit$surv      # survival curve for average patient
    ts <- fit$time
    nts <- length(ts)
    cmod <- coxph(formulaCen, x = TRUE, data = x)
    fit <- survfit(cmod)
    c0 <- fit$surv      # censoring curve for average patient

    # find event rate in base cohort (if everyone was followed to end of study)
    Xb <- as.vector(smod$x %*% coef(smod))
    mx <- colMeans(smod$x)
    xb0 <- mx %*% coef(smod)
    s0end <- min(s0)
    if(is.null(eventRate)) eventRate <- 1-mean(s0end^exp(Xb - xb0))

    # find delta value needed to get approximate desired event rate under new parameters
    bnew <- replace(MMOut*coef(smod), names(coef(smod)) == exposure, log(effectOR))
    Xbnew <- as.vector(smod$x %*% bnew)
    sXend <- s0end^(exp(Xb - xb0))
    fn <- function(d) mean(sXend^d) - (1 - eventRate)
    delta <- uniroot(fn, lower = 0, upper = 20)$root

    # setup n X nts matrix of individual survival and censoring curves under new parameters
    Sx <- matrix(unlist(lapply(s0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)
    Xbnew <- as.vector(smod$x %*% coef(cmod))
    xb0 <- mx %*% coef(cmod)
    Cx <- matrix(unlist(lapply(c0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)

    # fit PS model
    ## Estimate logit model for probability of exposure
    psmodel<- glm2(formulaExp, family = "binomial", data=x,control=glm.control(trace=TRUE))
    # find exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev <- mean(psmodel$y)
    ## Design matrix used for exposure logistic regression
    XEXP<- gam(formulaExp, data=x, family = "binomial", fit = FALSE)$X

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp<- c(coef(psmodel)[1], MMExp*coef(psmodel)[-1])
    XbnewExp<- as.vector(XEXP%*%bnewExp)
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    ids <- tnew <- ynew <-expnew<- data.frame(matrix(nrow = size, ncol = nsim))

    for(sim in 1:nsim) {
      idxs <- sample(n, size, replace = TRUE)
      ids[,sim] <- x$ID[idxs]

      # event time
      u <- runif(size, 0, 1)
      w <- apply(Sx[idxs,] < u, 1, function(x) which(x)[1]) # the first time survival drops below u
      stime <- ts[w]
      w <- Sx[idxs,nts] > u     # for any individuals with survival that never drops below u,
      stime[w] <- max(ts) + 1  # replace with arbitrary time beyond last observed event/censoring time

      # censoring time
      u <- runif(size, 0, 1)
      w <- apply(Cx[idxs,] < u, 1, function(x) which(x)[1]) # the first time censor-free survival drops below u
      ctime <- ts[w]
      w <- Cx[idxs,nts] > u     # for any individuals with censor-free survival that never drops below u,
      ctime[w] <- max(ts)    # replace with hard censor time at last observed event/censoring time


      expnew[,sim]<- rbinom(size,1,Probexp[idxs])

      # put it together
      tnew[,sim] <- pmin(stime, ctime)
      names(tnew) <- paste("TIME", 1:nsim, sep = "")
      ynew[,sim] <- stime == tnew[,sim]
      names(ynew) <- paste("EVENT", 1:nsim, sep = "")

    }

    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(tnew) <- paste("TIME", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,tnew,expnew)
    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp,Sim_Data = sim_out_bin))

  }

  else if(is.null(formulaOut)==FALSE & is.null(formulaCen)==FALSE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE & is.null(objectCen)==TRUE & is.null(objectExp)==TRUE){
    n <- nrow(data)
    exposure<-all.vars(formulaOut)[4]
    sidx <- sapply(c(idVar),function(v) which(names(data) == v))
    names(data)[sidx] <- c("ID")

    x<-data

    ##estimate survival and censoring models
    smod <- coxph(formulaOut, x = TRUE, data = x)
    fit <- survfit(smod)
    s0 <- fit$surv      # survival curve for average patient
    ts <- fit$time
    nts <- length(ts)
    cmod <- coxph(formulaCen, x = TRUE, data = x)
    fit <- survfit(cmod)
    c0 <- fit$surv      # censoring curve for average patient

    # find event rate in base cohort (if everyone was followed to end of study)
    Xb <- as.vector(smod$x %*% coef(smod))
    mx <- colMeans(smod$x)
    xb0 <- mx %*% coef(smod)
    s0end <- min(s0)
    if(is.null(eventRate)) eventRate <- 1-mean(s0end^exp(Xb - xb0))

    # find delta value needed to get approximate desired event rate under new parameters
    bnew <- replace(MMOut*coef(smod), names(coef(smod)) == exposure, log(effectOR))
    Xbnew <- as.vector(smod$x %*% bnew)
    sXend <- s0end^(exp(Xb - xb0))
    fn <- function(d) mean(sXend^d) - (1 - eventRate)
    delta <- uniroot(fn, lower = 0, upper = 20)$root

    # setup n X nts matrix of individual survival and censoring curves under new parameters
    Sx <- matrix(unlist(lapply(s0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)
    Xbnew <- as.vector(smod$x %*% coef(cmod))
    xb0 <- mx %*% coef(cmod)
    Cx <- matrix(unlist(lapply(c0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)


    #### sample and simulate
    ids <- tnew <- ynew <- data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs <- sample(n, size, replace = TRUE)
      ids[,sim] <- x$ID[idxs]

      # event time
      u <- runif(size, 0, 1)
      w <- apply(Sx[idxs,] < u, 1, function(x) which(x)[1]) # the first time survival drops below u
      stime <- ts[w]
      w <- Sx[idxs,nts] > u     # for any individuals with survival that never drops below u,
      stime[w] <- max(ts) + 1  # replace with arbitrary time beyond last observed event/censoring time

      # censoring time
      u <- runif(size, 0, 1)
      w <- apply(Cx[idxs,] < u, 1, function(x) which(x)[1]) # the first time censor-free survival drops below u
      ctime <- ts[w]
      w <- Cx[idxs,nts] > u     # for any individuals with censor-free survival that never drops below u,
      ctime[w] <- max(ts)    # replace with hard censor time at last observed event/censoring time


      # put it together
      tnew[,sim] <- pmin(stime, ctime)
      names(tnew) <- paste("TIME", 1:nsim, sep = "")
      ynew[,sim] <- stime == tnew[,sim]
      names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    }

    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(tnew) <- paste("TIME", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,tnew)
    return(list(TrueOutBeta = bnew, Sim_Data = sim_out_bin))
  }

  else if(is.null(formulaOut)==TRUE & is.null(formulaCen)==TRUE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE & is.null(objectCen)==TRUE & is.null(objectExp)==TRUE){
    n <- nrow(data)
    sidx <- sapply(c(idVar),function(v) which(names(data) == v))
    names(data)[sidx] <- c("ID")

    x<-data

    # fit PS model
    ## Estimate logit model for probability of exposure
    psmodel<- glm2(formulaExp, family = "binomial", data=x,control=glm.control(trace=TRUE))
    # find exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev <- mean(psmodel$y)
    ## Design matrix used for exposure logistic regression
    XEXP<- gam(formulaExp, data=x, family = "binomial", fit = FALSE)$X

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp<- c(coef(psmodel)[1], MMExp*coef(psmodel)[-1])
    XbnewExp<- as.vector(XEXP%*%bnewExp)
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    ids <- expnew<- data.frame(matrix(nrow = size, ncol = nsim))
    ##RR<-RD<- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      idxs <- sample(n, size, replace = TRUE)
      ids[,sim] <- x$ID[idxs]

      expnew[,sim]<- rbinom(size,1,Probexp[idxs])

    }
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids,expnew)
    return(list(TrueExpBeta = bnewExp,Sim_Data = sim_out_bin))

  }

  else if(is.null(formulaOut)==TRUE & is.null(formulaCen)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectCen)==FALSE & is.null(objectExp)==FALSE){

    dataOut<- as.data.frame(objectOut$x)
    n<- nrow(dataOut)
    names(dataOut)[1] <- "TREAT"
    OUTCOME<-as.numeric((objectOut$y)[,1])
    TIME<-as.numeric((objectOut$y)[,2])
    dataOut<- cbind(idVar,TIME,OUTCOME,dataOut)
    names(dataOut)[1]<- "ID"

    ##estimate survival and censoring models

    fit <- survfit(objectOut)
    s0 <- fit$surv      # survival curve for average patient
    ts <- fit$time
    nts <- length(ts)
    fit <- survfit(objectCen)
    c0 <- fit$surv      # censoring curve for average patient

    # find event rate in base cohort (if everyone was followed to end of study)
    Xb <- as.vector(objectOut$x %*% (coef(objectOut)))
    mx <- colMeans(objectOut$x)
    xb0 <- mx %*% coef(objectOut)
    s0end <- min(s0)
    if(is.null(eventRate)) eventRate <- 1-mean(s0end^exp(Xb - xb0))

    # find delta value needed to get approximate desired event rate under new parameters
    ##bnew <- replace(MMOut*coef(smod), names(coef(smod)) == exposure, log(effectOR))

    exposure<-names(coef(objectOut))[1]
    bnew <- replace(MMOut*coef(objectOut), names(coef(objectOut)) == exposure, log(effectOR))
    names(bnew)[1]<-"TREAT"
    Out<- as.data.frame(objectOut$x)
    names(Out)[1]<-"TREAT"
    Xbnew <- as.vector(as.matrix(Out)%*% bnew)
    sXend <- s0end^(exp(Xb - xb0))
    fn <- function(d) mean(sXend^d) - (1 - eventRate)
    delta <- uniroot(fn, lower = 0, upper = 20)$root

    # setup n X nts matrix of individual survival and censoring curves under new parameters
    Sx <- matrix(unlist(lapply(s0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)
    Xbnew <- as.vector(objectOut$x %*% coef(objectCen))
    xb0 <- mx %*% coef(objectCen)
    Cx <- matrix(unlist(lapply(c0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)

    # find exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev <- mean(objectExp$y)
    ## Design matrix used for exposure logistic regression
    XEXP<- model.matrix(objectExp)

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp<- c(coef(objectExp)[1], MMExp*coef(objectExp)[-1])
    XbnewExp<- as.vector(XEXP%*%bnewExp)
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    ids <- tnew <- ynew <-expnew<- data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs <- sample(n, size, replace = TRUE)
      ids[,sim] <- dataOut$ID[idxs]

      # event time
      u <- runif(size, 0, 1)
      w <- apply(Sx[idxs,] < u, 1, function(x) which(x)[1]) # the first time survival drops below u
      stime <- ts[w]
      w <- Sx[idxs,nts] > u     # for any individuals with survival that never drops below u,
      stime[w] <- max(ts) + 1  # replace with arbitrary time beyond last observed event/censoring time

      # censoring time
      u <- runif(size, 0, 1)
      w <- apply(Cx[idxs,] < u, 1, function(x) which(x)[1]) # the first time censor-free survival drops below u
      ctime <- ts[w]
      w <- Cx[idxs,nts] > u     # for any individuals with censor-free survival that never drops below u,
      ctime[w] <- max(ts)    # replace with hard censor time at last observed event/censoring time


      expnew[,sim]<- rbinom(size,1,Probexp[idxs])

      # put it together
      tnew[,sim] <- pmin(stime, ctime)
      names(tnew) <- paste("TIME", 1:nsim, sep = "")
      ynew[,sim] <- stime == tnew[,sim]
      names(ynew) <- paste("EVENT", 1:nsim, sep = "")

    }
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(tnew) <- paste("TIME", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,tnew,expnew)

    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp,Sim_Data = sim_out_bin))

  }

  else if(is.null(formulaOut)==TRUE & is.null(formulaCen)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectCen)==FALSE & is.null(objectExp)==TRUE){

    dataOut<- as.data.frame(objectOut$x)
    n<- nrow(dataOut)
    names(dataOut)[1] <- "TREAT"
    OUTCOME<-as.numeric((objectOut$y)[,1])
    TIME<-as.numeric((objectOut$y)[,2])
    dataOut<- cbind(idVar,TIME,OUTCOME,dataOut)
    names(dataOut)[1]<- "ID"

    ##estimate survival and censoring models

    fit <- survfit(objectOut)
    s0 <- fit$surv      # survival curve for average patient
    ts <- fit$time
    nts <- length(ts)
    fit <- survfit(objectCen)
    c0 <- fit$surv      # censoring curve for average patient

    # find event rate in base cohort (if everyone was followed to end of study)
    Xb <- as.vector(objectOut$x %*% (coef(objectOut)))
    mx <- colMeans(objectOut$x)
    xb0 <- mx %*% coef(objectOut)
    s0end <- min(s0)
    if(is.null(eventRate)) eventRate <- 1-mean(s0end^exp(Xb - xb0))

    # find delta value needed to get approximate desired event rate under new parameters
    ##bnew <- replace(MMOut*coef(smod), names(coef(smod)) == exposure, log(effectOR))

    exposure<-names(coef(objectOut))[1]
    bnew <- replace(MMOut*coef(objectOut), names(coef(objectOut)) == exposure, log(effectOR))
    names(bnew)[1]<-"TREAT"
    Out<- as.data.frame(objectOut$x)
    names(Out)[1]<-"TREAT"
    Xbnew <- as.vector(as.matrix(Out)%*% bnew)
    sXend <- s0end^(exp(Xb - xb0))
    fn <- function(d) mean(sXend^d) - (1 - eventRate)
    delta <- uniroot(fn, lower = 0, upper = 20)$root

    # setup n X nts matrix of individual survival and censoring curves under new parameters
    Sx <- matrix(unlist(lapply(s0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)
    Xbnew <- as.vector(objectOut$x %*% coef(objectCen))
    xb0 <- mx %*% coef(objectCen)
    Cx <- matrix(unlist(lapply(c0, function(s) s^(delta*exp(Xbnew - xb0)))), nrow = n)

    #### sample and simulate
    ids <- tnew <- ynew <- data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs <- sample(n, size, replace = TRUE)
      ids[,sim] <- dataOut$ID[idxs]

      # event time
      u <- runif(size, 0, 1)
      w <- apply(Sx[idxs,] < u, 1, function(x) which(x)[1]) # the first time survival drops below u
      stime <- ts[w]
      w <- Sx[idxs,nts] > u     # for any individuals with survival that never drops below u,
      stime[w] <- max(ts) + 1  # replace with arbitrary time beyond last observed event/censoring time

      # censoring time
      u <- runif(size, 0, 1)
      w <- apply(Cx[idxs,] < u, 1, function(x) which(x)[1]) # the first time censor-free survival drops below u
      ctime <- ts[w]
      w <- Cx[idxs,nts] > u     # for any individuals with censor-free survival that never drops below u,
      ctime[w] <- max(ts)    # replace with hard censor time at last observed event/censoring time


      # put it together
      tnew[,sim] <- pmin(stime, ctime)
      names(tnew) <- paste("TIME", 1:nsim, sep = "")
      ynew[,sim] <- stime == tnew[,sim]
      names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    }



    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(tnew) <- paste("TIME", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,tnew)
    return(list(TrueOutBeta = bnew,Sim_Data = sim_out_bin))



  }
  else if(is.null(formulaOut)==TRUE & is.null(formulaCen)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE & is.null(objectCen)==TRUE & is.null(objectExp)==FALSE){

    dataOut<- as.data.frame(model.matrix(objectExp)[,-1])
    n<- nrow(dataOut)
    TREAT<-objectExp$y
    dataOut<- cbind(idVar,TREAT,dataOut)
    names(dataOut)[1]<- "ID"


    # find exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev <- mean(objectExp$y)
    ## Design matrix used for exposure logistic regression
    XEXP<- model.matrix(objectExp)

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp<- c(coef(objectExp)[1], MMExp*coef(objectExp)[-1])
    XbnewExp<- as.vector(XEXP%*%bnewExp)
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    ids <-expnew<- data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs <- sample(n, size, replace = TRUE)
      ids[,sim] <- dataOut$ID[idxs]

      expnew[,sim]<- rbinom(size,1,Probexp[idxs])


    }
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids,expnew)

    return(list(TrueExpBeta = bnewExp,Sim_Data = sim_out_bin))


  }
}



