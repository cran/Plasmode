#' @title Performs the plasmode simulation
#' @description Creates 'plasmode'simulated datasets based on a given dataset when the outcome variable and exposure variable are binary. Plasmode simulation samples subjects with replacement from the observed data, uses subjects’ covariate data as is, and simulates exposure, outcome, or both.
#' @author Jessica M. Franklin, Younathan Abdia, and Shirley Wang
#' @param formulaOut An outcome model formula containing the binary outcome on the left-hand side and binary exposure along with potential confounders on the right-hand side.The functional form of the outcome model should be, Outcome ~ Exposure + Confounders. (Exposure main effect must be first independent variable)
#' @param objectOut A fitted model for the outcome model. The functional form of the fitted model for the outcome variable should be of form, Outcome ~ Exposure + Confounders.
#' @param formulaExp An exposure model formula containing the binary exposure on the left-hand side and potential confounders on the right-hand side. The functional form of the exposure model is, Exposure ~ Confounders.
#' @param objectExp A fitted model object for the exposure model.
#' @param data The dataset on which simulations are based. The data is required only when formulaOut or formulaExp or both are supplied to the argument.
#' @param idVar Name of the ID variable.
#' @param effectOR The desired treatment effect odds ratio. By default effectOR = 1.
#' @param MMOut A multiplier of confounder effects on outcome applied to the estimated log ORs in the outcome model. By default MMOut = 1 but one can specify a vector of length equivalent to the number of variables on the right-hand side of the outcome model.
#' @param MMExp A multiplier of confounder effects on exposure applied to the estimated log ORs in the exposure model. By default MMExp = 1 but one can specify a vector of length equivalent to the number of variables on the right-hand side of the exposure model.
#' @param nsim Number of desired simulated datasets.
#' @param size Desired size of simulated datasets (i.e., # of individuals).
#' @param eventRate Desired average event rate. Default is the event rate in the observed data.
#' @param exposedPrev Desired average exposure rate. Default is the esposure prevalence in the observed data.
#' @details At least one of formulaOut, formulaExp, objectOut, and objectExp must be specified, and which of these are specified will determine what gets simulated and how. If objectOut or objectExp are specified, these objects are used as the base model for outcome and exposure simulation. If formulaOut or formulaExp are specified, then data should be given and base models are fit in the data using glm2 with the given formulas. If formulaOut or objectOut is specified, outcome will be simulated based on subjects’ observed exposure. If formulaExp or objectExp is specified, exposure will be simulated. And if models are specified for both outcome and exposure, both variables will be simulated with simulated outcome dependent on the simulated exposure.
#' @export
#' @import mgcv
#' @import glm2
#' @import nlme
#' @importFrom stats as.formula coef fitted glm.control model.matrix plogis rbinom residuals rnorm runif uniroot var
#' @importFrom utils head tail
#' @return PlasmodeBin returns true beta coefficients used to generate the outcome and the exposure. It also returns the relative risk and risk difference estimated by the plasmode simulated data along with the data frame with the simulated data, including sampled IDs for each of nsim datasets along with simulated outcomes, exposure, or both.
#' \item{TrueOutBeta}{True beta coefficients used to generate the outcome.}
#' \item{TrueExpBeta}{True beta coefficients used to generate the exposure.}
#' \item{RR}{True relative risk estimated using the plasmode simulated data.}
#' \item{RD}{True risk difference estimated using the plasmode simulated data.}
#' \item{Sim_Data}{Plasmode simulated data, including sampled IDs for each of nsim datasets along with simulated outcomes, exposure, or both.}
#' @examples {
#' library(mgcv)
#' library(nlme)
#' library(glm2)
#' library(arm)
#' library(MASS)
#' library(lme4)
#' library(epiDisplay)
#' library(foreign)
#' library(nnet)
#'
#' data("Compaq")
#' levels(Compaq$stage) <- c(1,2,3,4)
#' Compaq$stage<-as.numeric(levels(Compaq$stage))[Compaq$stage]
#' ## Creating the binary exposure variable
#' Compaq$exposure<-ifelse(Compaq$hospital == "Public hospital",1,0)
#' ## Creating binary variables for some confounders
#' Compaq$ses1<-ifelse(Compaq$ses == "Poor",1,0)
#' Compaq$ses2<-ifelse(Compaq$ses == "Poor-middle",1,0)
#' Compaq$ses3<-ifelse(Compaq$ses == "High-middle",1,0)
#'
#' Compaq$age1<-ifelse(Compaq$agegr == "<40",1,0)
#' Compaq$age2<-ifelse(Compaq$agegr == "40-49",1,0)
#' Compaq$age3<-ifelse(Compaq$agegr == "50-59",1,0)
#'
#'
#' ## Creating the formulas for the outcome and the exposure model
#' form1<- status~ exposure + stage + ses1 + ses2 + ses3 + age1 + age2 + age3
#' form2<- exposure ~ stage + ses1 + ses2 + ses3 + age1 + age2 + age3
#'
#' set.seed(111)
#'Bin_Form1<-PlasmodeBin(formulaOut=form1, objectOut=NULL,formulaExp=form2,
#'                      objectExp= NULL,data=Compaq,idVar="id",effectOR =1,
#'                      MMOut=c(1,1,2,1,1,2,1,2),MMExp=c(1,1,1,1,1,1,1),
#'                      nsim=2, size=nrow(Compaq), eventRate=NULL, exposedPrev=NULL)
#'
#' Bin_Form2<-PlasmodeBin(formulaOut=form1, objectOut=NULL,formulaExp=NULL,
#'                       objectExp= NULL,data=Compaq,idVar="id",effectOR =1,
#'                       MMOut=c(1,1,2,1,1,2,1,2),MMExp=1, nsim=2,
#'                       size=nrow(Compaq), eventRate=NULL, exposedPrev=NULL)
#'
#' Bin_Form3<-PlasmodeBin(formulaOut=NULL, objectOut=NULL,formulaExp=form2,
#'                       objectExp= NULL,data=Compaq,idVar="id",effectOR =1,
#'                       MMOut=1,MMExp=c(1,1,1,1,1,1,1), nsim=2,
#'                       size=nrow(Compaq), eventRate=NULL, exposedPrev=NULL)
#'
#'
#' ###################################################################################################
#' ## One can provide the fitted model for the outcome model and the exposure model estimated by
#' ## glm, gam, and bayesglm. The functional form of the fitted model for the outcome variable should
#' ## of the form Outcome ~ Exposure + Confounders. The functional form of the exposure model is,
#' ## Exposure ~ Confounders.
#' ####################################################################################################
#'
#' Coeff1<- bayesglm(form1, family = "binomial", data=Compaq,control=glm.control(trace=TRUE))
#' Coeff2<- bayesglm(form2, family = "binomial", data=Compaq,control=glm.control(trace=TRUE))
#' sizesim<-nrow(model.matrix(Coeff1))
#' sizesim1<-nrow(model.matrix(Coeff2))
#'
#' Bin_Obj1<-PlasmodeBin(formulaOut=NULL, objectOut=Coeff1,formulaExp=NULL,
#'                      objectExp = Coeff2, idVar=Compaq$id,effectOR =1,
#'                      MMOut=c(1.5,1,2,1,1,1,1,1),MMExp=c(1,1,1,1,1,1,1),
#'                      nsim=2, size=sizesim, eventRate=NULL, exposedPrev=NULL)
#'
#' Bin_Obj2<-PlasmodeBin(formulaOut=NULL, objectOut=Coeff1,formulaExp=NULL,
#'                      objectExp = NULL,idVar=Compaq$id,effectOR =1,
#'                      MMOut=c(1.5,1,2,1,1,1,1,1),MMExp=1,
#'                      nsim=2, size=sizesim, eventRate=NULL, exposedPrev=NULL)
#'
#' Bin_Obj3<-PlasmodeBin(formulaOut=NULL, objectOut=NULL,formulaExp=NULL,
#'                      objectExp = Coeff2,idVar=Compaq$id,effectOR =1, MMOut=1,
#'                      MMExp=c(1,1,1,1,1,1,1),
#'                      nsim=2, size=sizesim1, eventRate=NULL, exposedPrev=NULL)
#' }




PlasmodeBin<- function(formulaOut=NULL, objectOut=NULL,formulaExp=NULL,objectExp=NULL,data, idVar,
                       effectOR =1, MMOut=1,MMExp=1, nsim, size, eventRate=NULL, exposedPrev=NULL)

{

  ## Code for simulating data when the outcome is binary and data set is provided to estimated the outcome and exposure.
  if(is.null(formulaOut)==FALSE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE & is.null(objectExp)==TRUE)

  {
    outcome<- all.vars(formulaOut)[1] ## selects the outcome variable
    exposure<- all.vars(formulaOut)[2] ##selects the exposure variable

    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # estimate logit model for probability of outcome
    modOutBin<- glm2(formulaOut, family = "binomial", data=x,control=glm.control(trace=TRUE))
    ## Design Matrix used for outcome logistic regression
    X <- gam(formulaOut, data=x, family = "binomial", fit = FALSE)$X
    pred <- fitted(modOutBin) # linear predictor


    # find event rate in base cohort
    if(is.null(eventRate)) eventRate <- mean(pred)

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(coef(modOutBin)[1], MMOut*coef(modOutBin)[-1])
    bnew <- replace(bnew, names(coef(modOutBin)) == exposure, log(effectOR))
    Xbnew <- as.vector(X %*% bnew)
    fn <- function(d) mean(1 - 1/(1 + exp(d+Xbnew))) - eventRate
    delta <- uniroot(fn, lower = -20, upper = 20)$root
    pnew <- 1 - 1/(1 + exp(delta+Xbnew))


    ## Estimate logit model for probability of exposure
    modExp<- glm2(formulaExp, family = "binomial", data,control=glm.control(trace=TRUE))
    ## Design matrix used for exposure logistic regression
    XEXP<- gam(formulaExp, data, family = "binomial", fit = FALSE)$X

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev))exposedPrev<- mean(modExp$y)

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp<- c(coef(modExp)[1], MMExp*coef(modExp)[-1])
    XbnewExp<- as.vector(XEXP%*%bnewExp)
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)


    #### sample and simulate
    ids <- ynew <- expnew<-data.frame(matrix(nrow = size, ncol = nsim))
    RR<-RD<- vector('numeric', length = nsim)
    ##datasim<-list()
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      ynew[,sim] <- rbinom(size, 1, pnew[c(idxs0,idxs1)])
      expnew[,sim]<- rbinom(size,1,Probexp[c(idxs0,idxs1)])
      datasim<-X[c(idxs0,idxs1),]
      datasim[,2]<-1
      p_1<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      datasim[,2]<-0
      p_0<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)


    ## Finding the simulated treatment at new exposure prevalence.

    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and data set is provided to estimated the outcome.

  else if(is.null(formulaOut)==FALSE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE & is.null(objectExp)==TRUE)

  {
    outcome<- all.vars(formulaOut)[1] ## selects the outcome variable
    exposure<- all.vars(formulaOut)[2] ##selects the exposure variable

    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # estimate logit model for probability of outcome
    modOutBin<- glm2(formulaOut, family = "binomial", data=x,control=glm.control(trace=TRUE))
    ## Design Matrix used for outcome logistic regression
    X <- gam(formulaOut, data=x, family = "binomial", fit = FALSE)$X
    pred <- fitted(modOutBin) # linear predictor

    # find event rate in base cohort
    if(is.null(eventRate)) eventRate <- mean(pred)

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(coef(modOutBin)[1], MMOut*coef(modOutBin)[-1])
    bnew <- replace(bnew, names(coef(modOutBin)) == exposure, log(effectOR))
    Xbnew <- as.vector(X %*% bnew)
    fn <- function(d) mean(1 - 1/(1 + exp(d+Xbnew))) - eventRate
    delta <- uniroot(fn, lower = -20, upper = 20)$root
    pnew <- 1 - 1/(1 + exp(delta+Xbnew))
    rm(modOutBin)

    ### sample and simulate
    ids <- ynew <-data.frame(matrix(nrow = size, ncol = nsim))
    RR<-RD<- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      ynew[,sim] <- rbinom(size, 1, pnew[c(idxs0,idxs1)])
      datasim<-X[c(idxs0,idxs1),]
      datasim[,2]<-1
      p_1<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      datasim[,2]<-0
      p_0<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }

    ARR<-mean(RR)
    ARD<-mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew)


    ## Finding the simulated treatment at new exposure prevalence.

    return(list(TrueOutBeta = bnew, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and data set is provided to estimated the exposure.
  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE & is.null(objectExp)==TRUE)

  {
    exposure<- all.vars(formulaExp)[1] ##selects the exposure variable

    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    ## Estimate logit model for probability of exposure
    modExp<- glm2(formulaExp, family = "binomial", data,control=glm.control(trace=TRUE))
    ## Design matrix used for exposure logistic regression
    XEXP<- gam(formulaExp, data, family = "binomial", fit = FALSE)$X

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev))exposedPrev<- mean(modExp$y)
    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp<- c(coef(modExp)[1], MMExp*coef(modExp)[-1])
    XbnewExp<- as.vector(XEXP%*%bnewExp)
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)
    rm(modExp, XEXP)

    #### sample and simulate
    ids <- expnew<-data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      expnew[,sim]<- rbinom(size,1,Probexp[c(idxs0,idxs1)])
    }

    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    ##save(ids,ynew,file="simdata.RData")
    sim_out_bin<-data.frame(ids,expnew)


    ## Finding the simulated treatment at new exposure prevalence.

    return(list(TrueExpBeta = bnewExp, Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and objects are provided to estimated the outcome and exposure.

  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectExp)==FALSE)

  {
    DesMatOut<- model.matrix(objectOut)
    exposure<- colnames(DesMatOut)[2] ##selects the exposure variable
    dataOut<- as.data.frame(DesMatOut)
    dataOut<- cbind(idVar,dataOut)

    x <- dataOut[order(dataOut[,exposure]),]# order according to exposure status, unexposed first
    X<-x[-1]

    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # find event rate in base cohort
    OutCoeff<-coef(objectOut)
    if(is.null(eventRate)) eventRate <- mean(1-1/(1+exp(OutCoeff%*%t(X))))

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(OutCoeff[1], MMOut*OutCoeff[-1])
    bnew <- replace(bnew, names(OutCoeff) == exposure, log(effectOR))
    Xbnew <- as.vector(bnew%*%t(X))
    fn <- function(d) mean(1 - 1/(1 + exp(d+Xbnew))) - eventRate
    delta <- uniroot(fn, lower = -20, upper = 20)$root
    pnew <- 1 - 1/(1 + exp(delta+Xbnew))

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev))exposedPrev<- sum(X[2])/nrow(X)
    DesMatExp<-model.matrix(objectExp)
    XEXP<- as.data.frame(DesMatExp)

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    ExpCoeff<-coef(objectExp)
    bnewExp<- c(ExpCoeff[1], MMExp*ExpCoeff[-1])
    XbnewExp<- as.vector(bnewExp%*%t(XEXP))
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    idvar<-colnames(dataOut)[1]
    ids <- ynew <- expnew<-data.frame(matrix(nrow = size, ncol = nsim))
    RR<-RD<- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- dataOut[idxs0, idvar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- dataOut[idxs1, idvar]
      ynew[,sim] <- rbinom(size, 1, pnew[c(idxs0,idxs1)])
      expnew[,sim]<- rbinom(size,1,Probexp[c(idxs0,idxs1)])
      datasim<-DesMatOut[c(idxs0,idxs1),]
      datasim[,2]<-1
      p_1<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      datasim[,2]<-0
      p_0<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)
    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectExp)==TRUE)

  {
    DesMatOut<- model.matrix(objectOut)
    exposure<- colnames(DesMatOut)[2] ##selects the exposure variable
    dataOut<- as.data.frame(DesMatOut)
    dataOut<- cbind(idVar,dataOut)

    x <- dataOut[order(dataOut[,exposure]),]# order according to exposure status, unexposed first
    X<-x[-1]

    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # find event rate in base cohort
    OutCoeff<-coef(objectOut)
    if(is.null(eventRate)) eventRate <- mean(1-1/(1+exp(OutCoeff%*%t(X))))

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(OutCoeff[1], MMOut*OutCoeff[-1])
    bnew <- replace(bnew, names(OutCoeff) == exposure, log(effectOR))
    Xbnew <- as.vector(bnew%*%t(X))
    fn <- function(d) mean(1 - 1/(1 + exp(d+Xbnew))) - eventRate
    delta <- uniroot(fn, lower = -20, upper = 20)$root
    pnew <- 1 - 1/(1 + exp(delta+Xbnew))

    #### sample and simulate
    idvar<-colnames(dataOut)[1]
    ids <- ynew <- data.frame(matrix(nrow = size, ncol = nsim))
    RR<-RD<- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- dataOut[idxs0, idvar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- dataOut[idxs1, idvar]
      ynew[,sim] <- rbinom(size, 1, pnew[c(idxs0,idxs1)])
      datasim<-DesMatOut[c(idxs0,idxs1),]
      datasim[,2]<-1
      p_1<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      datasim[,2]<-0
      p_0<- 1 - 1/(1+exp(as.vector(datasim %*% bnew)))
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew)
    return(list(TrueOutBeta = bnew, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and objects are provided to estimated the exposure.

  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE & is.null(objectExp)==FALSE)

  {
    DesMatExp<-model.matrix(objectExp)
    dataOut<- as.data.frame(DesMatExp)
    ex<- as.data.frame(objectExp$y)
    names(ex)<- "expo"
    dataOut<- cbind(idVar,ex,dataOut)
    exposure<- names(dataOut)[2]
    x <- dataOut[order(dataOut[,exposure]),][-1] # order according to exposure status, unexposed first

    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev))exposedPrev<- sum(x[1])/nrow(x)
    XEXP<- as.data.frame(DesMatExp)

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    ExpCoeff<-coef(objectExp)
    bnewExp<- c(ExpCoeff[1], MMExp*ExpCoeff[-1])
    XbnewExp<- as.vector(bnewExp%*%t(XEXP))
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    idvar<-names(dataOut)[1]
    ids <- expnew<-data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- dataOut[idxs0, idvar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- dataOut[idxs1, idvar]
      expnew[,sim]<- rbinom(size,1,Probexp[c(idxs0,idxs1)])
    }

    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, expnew)
    return(list(TrueExpBeta = bnewExp, Sim_Data = sim_out_bin))
  }

}



