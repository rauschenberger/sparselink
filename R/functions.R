
if(FALSE){
  rm(list=ls(all.names=TRUE))
  #install.packages(c("roxygen2","pkgdown","rcmdcheck","usethis","remotes","testthat","devtools"))
  setwd("C:/Users/arauschenberger/Desktop/sparselink/package")
  roxygen2::roxygenise()
  # usethis::use_mit_license()
  rcmdcheck::rcmdcheck()
  pkgdown::check_pkgdown()
  devtools::build()
  devtools::submit_cran()
  #pkgdown::build_site()
}

#----- main functions -----

#'@title
#'Sparse regression for related problems
#'
#'@description
#'Estimates sparse regression models (i.e., performing feature selection)
#'in multi-task learning or transfer learning.
#'Multi-task learning involves multiple targets,
#'and transfer learning involves multiple datasets.
#'
#'@export
#'
#'@param x \eqn{n \times p} matrix (multi-task learning) 
#'or list of \eqn{n_k \times p} matrices (transfer learning)
#'@param y \eqn{n \times q} matrix (multi-task learning)
#'or list of \eqn{n_k}-dimensional vectors (transfer learning)
#'@param family character \code{"gaussian"} or \code{"binomial"}
#'@param alpha.init elastic net mixing parameter for initial regressions,
#'default: 0.95 (lasso-like elastic net)
#'@param alpha elastic net mixing parameter of final regressions,
#'default: 1 (lasso)
#'@param nfolds number of internal cross-validation folds,
#'default: 10 (10-fold cross-validation)
#'@param type default \code{"exp"} scales weights with
#'\eqn{w_{ext}^{v_{ext}}+w_{int}^{v_{int}}}
#'(see internal function \link{construct_penfacs} for details)
#'@param cands candidate values for both scaling parameters,
#'default: \code{NULL} (\{0, 0.2, 0.4, 0.6, 0.8, 1\})
#'
#'@return
#'Returns an object of class \code{sparselink}, a list with multiple slots:
#'\itemize{
#'\item Stage 1 regressions (before sharing information):
#'Slot \code{glm.one} contains \eqn{q} objects of type
#'\code{cv.glmnet} (one for each problem).
#'\item Candidate scaling parameters (exponents):
#'Slot \code{weight} contains
#'a data frame with \eqn{n} combinations of exponents for the external (source)
#'and internal (target) weights
#'\item Stage 2 regressions (after sharing information):
#'Slot \code{glm.two} contains \eqn{q} lists (one for each problem)
#'of \eqn{n} objects of type \code{cv.glmnet}
#'(one for each combination of exponents).
#'\item Optimal regularisation parameters: Slot \code{lambda.min}
#'contains the cross-validated regularisation parameters for the stage 2
#'regressions.
#'\item Optimal scaling parameters: Slots \code{weight.ind} and
#'\code{weight.min} indicate or contain the cross-validated scaling parameters.
#'}
#'
#'@inherit sparselink-package references
#'
#'@seealso
#'Use \code{\link[=coef.sparselink]{coef}} to extract coefficients
#'and \code{\link[=predict.sparselink]{predict}} to make predictions.
#'
#'@examples
#'#--- multi-task learning ---
#'n <- 100
#'p <- 200
#'q <- 3
#'\dontshow{n <- 10; p <- 5; q <- 2}
#'family <- "gaussian"
#'x <- matrix(data=rnorm(n=n*p),nrow=n,ncol=p)
#'y <- matrix(data=rnorm(n*q),nrow=n,ncol=q)
#'object <- sparselink(x=x,y=y,family=family)
#'
#'#--- transfer learning ---
#'n <- c(100,50)
#'p <- 200
#'\dontshow{n <- c(10,10); p <- 5}
#'x <- lapply(X=n,function(x) matrix(data=stats::rnorm(n*p),nrow=x,ncol=p))
#'y <- lapply(X=n,function(x) stats::rnorm(x))
#'family <- "gaussian"
#'object <- sparselink(x=x,y=y,family=family)
#'
sparselink <- function(x,y,family,alpha.init=0.95,alpha=1,type="exp",nfolds=10,cands=NULL){
  
  alpha.one <- alpha.init
  alpha.two <- alpha
  
  cat(paste0("alpha.init=",alpha.init,", alpha=",alpha,", type=",type,"\n"))
  
  if(is.matrix(y) & is.matrix(x)){
    message("mode: multi-target learning")
    p <- ncol(x)
    q <- ncol(y)
    n <- rep(x=nrow(x),times=q)
    foldid <- make_folds_multi(y=y,family=family,nfolds=nfolds)
    y <- apply(y,2,function(x) x,simplify=FALSE)
    x <- replicate(n=q,expr=x,simplify=FALSE)
    foldid <- replicate(n=q,expr=foldid,simplify=FALSE)
    mode <- "multiple"
  } else if(is.list(y) & is.list(x)){
    message("mode: transfer learning")
    n <- sapply(X=y,FUN=base::length)
    p <- ncol(x[[1]])
    q <- length(x)
    foldid <- make_folds_trans(y=y,family=family,nfolds=nfolds)
    mode <- "transfer"
  } else {
    stop("Provide both x and y either as matrices (multi-target learning) or lists (transfer learning).")
  }
  
  if(length(family)==1){
    family <- rep(family,times=q)
  }
  
  # standardisation (for initial regressions)
  y_sta <- x_sta <- list()
  for(i in seq_len(q)){
    if(family[i]=="gaussian"){
      y_sta[[i]] <- scale(y[[i]],center=TRUE,scale=TRUE)
    } else {
      y_sta[[i]] <- y[[i]]
    }
    x_sta[[i]] <- scale(x[[i]],center=TRUE,scale=TRUE)
    x_sta[[i]][is.na(x_sta[[i]])] <- 0
  }
  
  # overall fit
  glm.one.ext <- glm.two.ext <- list()
  beta.ext <- matrix(data=NA,nrow=p,ncol=q)
  for(i in seq_len(q)){
    if(is.na(alpha.one)){
      beta.ext[,i] <- stats::cor(x=x[[i]],y=y[[i]],method="spearman")
      beta.ext[,i][is.na(beta.ext[,i])] <- 0
    } else {
      glm.one.ext[[i]] <- glmnet::cv.glmnet(x=x_sta[[i]],y=y_sta[[i]],family=family[i],alpha=alpha.one,foldid=foldid[[i]],standardize=FALSE)
      beta.ext[,i] <- stats::coef(object=glm.one.ext[[i]],s="lambda.min")[-1]
    }
  }
  
  if(type %in% c("geo","exp","rem","ari")){
    if(is.null(cands)){
      cands <- seq(from=0,to=1,by=0.2)
    }
    weight <- expand.grid(source=cands,target=cands)
  } else if(type %in% c("geo.con","exp.con","rem.con","ari.con")){
    cands <- seq(from=0,to=1,by=0.05)
    weight <- data.frame(source=cands,target=1-cands)
  } else {
    stop(paste0("Invalid type=",type))
  }
  
  for(i in seq_len(q)){
    rest.ext <- construct_weights(coef=beta.ext,id=i)
    glm.two.ext[[i]] <- list()
    for(l in seq_len(nrow(weight))){
      pf <- construct_penfacs(w_int=rest.ext$weight.target,w_ext=rest.ext$weight.source,v_int=weight$target[l],v_ext=weight$source[l],type=type)
      if(all(pf==Inf)){
        glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family[i],lambda=99e99,alpha=alpha.two)
      } else {
        glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family[i],lower.limits=rest.ext$lower.limits,upper.limits=rest.ext$upper.limits,penalty.factor=pf,alpha=alpha.two)
      }
    }
  }
  
  # cross-validation
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- list()
    for(l in seq_len(nrow(weight))){
      nlambda <- length(glm.two.ext[[i]][[l]]$lambda)
      y_hat[[i]][[l]] <- matrix(data=NA,nrow=n[i],ncol=nlambda)
    }
  }
  
  for(k in seq_len(nfolds)){
    if(mode=="multiple"){
      beta.int <- matrix(data=NA,nrow=p,ncol=q) 
      for(i in seq_len(q)){
        cond <- foldid[[i]]==k # switch to vector?
        if(is.na(alpha.one)){
          beta.int[,i] <- stats::cor(x=x[[i]][!cond,],y=y[[i]][!cond],method="spearman")
          beta.int[,i][is.na(beta.int[,i])] <- 0
        } else {
          glm.one.int <- glmnet::glmnet(x=x_sta[[i]][!cond,],y=y_sta[[i]][!cond],family=family[i],alpha=alpha.one,standardize=FALSE)
          beta.int[,i] <- stats::coef(object=glm.one.int,s=glm.one.ext[[i]]$lambda.min)[-1]
        }
      }
    }
    for(i in seq_len(q)){
      if(mode=="transfer"){
        beta.int <- beta.ext
        cond <- foldid[[i]]==k
        if(is.na(alpha.one)){
          beta.int[,i] <- stats::cor(x=x[[i]][!cond,],y=y[[i]][!cond],method="spearman")
          beta.int[,i][is.na(beta.int[,i])] <- 0
        } else {
          glm.one.int <- glmnet::glmnet(x=x_sta[[i]][!cond,],y=y_sta[[i]][!cond],family=family[i],alpha=alpha.one,standardize=FALSE)
          beta.int[,i] <- stats::coef(object=glm.one.int,s=glm.one.ext[[i]]$lambda.min)[-1]
        }
      }
      rest.int <- construct_weights(coef=beta.int,id=i)
      for(l in seq_len(nrow(weight))){
        pf <- construct_penfacs(w_int=rest.int$weight.target,w_ext=rest.int$weight.source,v_int=weight$target[l],v_ext=weight$source[l],type=type)
        if(all(pf==Inf)){
          glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family[i],lambda=99e99,alpha=alpha.two)
        } else {
          glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family[i],lower.limits=rest.int$lower.limits,upper.limits=rest.int$upper.limits,penalty.factor=pf,alpha=alpha.two)
        }
        y_hat[[i]][[l]][cond,] <- stats::predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=glm.two.ext[[i]][[l]]$lambda,type="response")
      }
    }
  }
  
  metric <- list()
  cvm.min <- lambda.ind <- lambda.min <- matrix(data=NA,nrow=nrow(weight),ncol=q)
  for(i in seq_len(q)){
    metric[[i]] <- list()
    for(l in seq_len(nrow(weight))){
      metric[[i]][[l]] <- apply(X=y_hat[[i]][[l]],MARGIN=2,FUN=function(x)
        calc_metric(y=y[[i]],y_hat=x,family=family[i]))
      lambda.ind[l,i] <- which.min(metric[[i]][[l]])
      lambda.min[l,i] <- glm.two.ext[[i]][[l]]$lambda[lambda.ind[l,i]]
    }
    cvm.min[,i] <- sapply(X=metric[[i]],FUN=min)
    #if(trial){
    #  tryCatch(expr=plot_weight(x=weight,y=cvm.min[,i]),error=function(x) NULL)
    #}
  }
  weight.ind <- apply(X=cvm.min,MARGIN=2,FUN=which.min)
  weight.min <- weight[weight.ind,]
  
  list <- list(glm.one=glm.one.ext,glm.two=glm.two.ext,weight=weight,weight.ind=weight.ind,weight.min=weight.min,lambda.min=lambda.min,info=data.frame(p=p,q=q,mode=mode,family=paste0(family,collapse=", ")))
  class(list) <- "sparselink"
  return(list)
}

#'@title Print \code{sparselink} object
#'@description
#'Prints an object of class \code{sparselink}
#'
#'@export
#'@keywords internal
#'
#'@param x object of class \code{sparselink}
#'(generated by function \link{sparselink})
#'@inheritParams predict.sparselink
#'
#'@examples
#'n <- 100; p <- 50; q <- 3
#'family <- "gaussian"
#'x <- matrix(data=rnorm(n=n*p),nrow=n,ncol=p)
#'y <- matrix(data=rnorm(n*q),nrow=n,ncol=q)
#'object <- sparselink(x=x,y=y,family=family)
#'object
#'
print.sparselink <- function(x,...){
  if(x$info$mode=="multiple"){
    mode <- "MTL"
  } else if(x$info$mode=="transfer"){
    mode <- "TL"
  }
  cat("sparselink-object:\n")
  cat(ifelse(mode=="MTL","multi-task","transfer"),"learning\n")
  cat(paste0(x$info$q," problems (",ifelse(mode=="MTL","targets","datasets"),")\n"))
  cat(x$info$p,"features")
}

#'@title
#'Regression Coefficients
#'
#'@description
#'Extracts coefficients from multi-task or transfer learning regression model.
#'
#'@export
#'
#'@inheritParams predict.sparselink
#'
#'@return
#'Returns estimated coefficients.
#'The output is a list with two slots:
#'slot \code{alpha} with the estimated intercept
#'(vector of length \eqn{q}),
#'and slot \code{beta} with the estimated slopes
#'(matrix with \eqn{p} rows and \eqn{q} columns).
#'
#'@inherit sparselink-package references
#'
#'@seealso
#'Use \code{\link{sparselink}} to fit the model
#'and \code{\link[=predict.sparselink]{predict}} to make predictions.
#'
#'@examples
#'family <- "gaussian"
#'data <- sim_data_trans(family=family)
#'\dontshow{data <- sim_data_trans(family=family,n0=10,p=3)}
#'#data <- sim_data_multi(family=family)
#'object <- sparselink(x=data$X_train,y=data$y_train,family=family)
#'coef <- coef(object=object)
#'
coef.sparselink <- function(object,...){
  id <- object$weight.ind
  p <- object$info$p
  q <- object$info$q
  alpha <- rep(x=NA,times=object$info$q)
  beta <- matrix(data=NA,nrow=object$info$p,ncol=object$info$q)
  for(i in seq_len(q)){
    temp <- stats::coef(object=object$glm.two[[i]][[id[i]]],s=object$lambda.min[id[i],i])
    alpha[i] <- temp[1]
    beta[,i] <- temp[-1][seq_len(p)]+temp[-1][seq(from=p+1,to=2*p)]
  }
  coef <- list(alpha=alpha,beta=beta)
  return(coef)
}

#'@title
#'Out-of-sample Predictions
#'
#'@description
#'Predicts outcomes with a multi-task or transfer learning regression model.
#'
#'@export
#'
#'@param object
#'object of class \code{"sparselink"}
#'(generated by function \link{sparselink})
#'
#'@param newx
#'features:
#'matrix with \eqn{n} rows (samples) and \eqn{p} columns (variables)
#'for multi-task learning;
#'list of \eqn{q} matrices
#'with \eqn{n_k} rows (samples) and \eqn{p} columns (variables)
#'for transfer learning, for each \eqn{k} in \eqn{1,\ldots,q}
#'
#'@param weight
#'hyperparameters for scaling external and internal weights:
#'numeric vector of length 2,
#'with the first entry for the external weights
#'(prior coefficients from source data),
#'and the second entry for the internal weights
#'(prior coefficients from target data),
#'selected values must be among the candidate values,
#'default: \code{NULL} (using cross-validated weights)
#'
#'@param ... (not applicable)
#'
#'@return
#'Returns predicted values or predicted probabilities.
#'The output is a list of \eqn{q} column vectors of length \eqn{n_k}
#'for \eqn{k} in \eqn{1,\ldots,q}.
#'Each vector corresponds to one target (multi-task learning)
#'or one dataset (transfer learning).
#'
#'@inherit sparselink-package references
#'
#'@seealso
#'Use \code{\link{sparselink}} to fit the model
#'and \code{\link[=coef.sparselink]{coef}} to extract coefficients.
#'
#'@examples
#'family <- "gaussian"
#'data <- sim_data_multi(family=family)
#'\dontshow{data <- sim_data_multi(family=family,n0=10,p=3)}
#'#data <- sim_data_trans(family=family)
#'object <- sparselink(x=data$X_train,y=data$y_train,family=family)
#'y_hat <- predict(object=object,newx=data$X_test)
#'
predict.sparselink <- function(object,newx,weight=NULL,...){
  if(is.null(weight)){
    id <- object$weight.ind
  } else {
    if(length(weight)!=length(object$glm.one)){stop("invalid length")}
    id <- which(object$weight$source==weight[1] & object$weight$target==weight[2])
  }
  y_hat <- list()
  if(is.list(newx)){
    q <- length(newx)
    for(i in seq_len(q)){
      y_hat[[i]] <- stats::predict(object=object$glm.two[[i]][[id[i]]],newx=cbind(newx[[i]],newx[[i]]),s=object$lambda.min[id[i],i],type="response")
    }
  } else {
    for(i in seq_len(object$info$q)){
      y_hat[[i]] <- stats::predict(object=object$glm.two[[i]][[id[i]]],newx=cbind(newx,newx),s=object$lambda.min[id[i],i],type="response")
    }
  }
  return(y_hat)
}

#----- internal functions -----

#'@title logit function
#'@export
#'@keywords internal
#'
#'@param x numeric vector with values in unit interval
#'
#'@examples
#'x <- seq(from=0,to=1,length.out=100)
#'y <- logit(x=x)
#'graphics::plot(x=x,y=y,type="l")
#'graphics::abline(v=0.5,lty=2)
#'graphics::abline(h=0,lty=2)
#'
logit <- function(x){
  if(any(x<0|x>1)){stop("x must be in unit interval")}
  log(x/(1-x))
}

#'@title Sigmoid function
#'@export
#'@keywords internal
#'
#'@param x numeric vector
#'
#'@examples
#'x <- seq(from=-3,to=3,length.out=100)
#'y <- sigmoid(x)
#'graphics::plot(x=x,y=y,type="l")
#'graphics::abline(v=0,lty=2)
#'graphics::abline(h=0.5,lty=2)
#'
sigmoid <- function(x){
  1/(1+exp(-x))
}

#'@title Link function
#'@export
#'@keywords internal
#'
#'@description
#'Applies the link function.
#'
#'@param mu numeric vector (with values in unit interval if family="binomial")
#'@param family character \code{"gaussian"} or \code{"binomial"}
#'
#'@examples
#'family <- "binomial"
#'from <- ifelse(family=="binomial",0,-3)
#'to <- ifelse(family=="binomial",1,3)
#'mu <- seq(from=from,to=to,length.out=100)
#'eta <- link_function(mu=mu,family=family)
#'graphics::plot(x=mu,y=eta,type="l",main=family)
#'v <- ifelse(family=="binomial",0.5,0)
#'graphics::abline(v=v,lty=2)
#'graphics::abline(h=0,lty=2)
#'
link_function <- function(mu,family){
  if(family=="gaussian"){
    eta <- mu
  } else if(family=="binomial"){
    eta <- logit(mu)
  } else {
    stop("family=\"",family,"\" not implemented")
  }
  return(eta)
}

#'@title Mean function
#'@export
#'@keywords internal
#'
#'@description
#'Applies the mean function (inverse link function).
#'
#'@param eta numeric vector
#'@param family character \code{"gaussian"} or \code{"binomial"}
#'
#'@examples
#'family <- "binomial"
#'eta <- seq(from=-3,to=3,length.out=100)
#'mu <- mean_function(eta=eta,family=family)
#'graphics::plot(x=eta,y=mu,type="l",main=family)
#'graphics::abline(v=0,lty=2)
#'h <- ifelse(family=="binomial",0.5,0)
#'graphics::abline(h=h,lty=2)
#'
mean_function <- function(eta,family){
  if(family=="gaussian"){
    mu <- eta
  } else if(family=="binomial"){
    mu <- sigmoid(eta)
  } else {
    stop("family=\"",family,"\" not implemented")
  }
  return(mu)
}

#'@title Calculate deviance
#'@export
#'@keywords internal
#'
#'@description
#'Calculates Gaussian deviance (mean-squared error) and binomial deviance.
#'
#'@param y response
#'@param y_hat predictor
#'@param family character \code{"gaussian"} or \code{"binomial"}
#'
#'@examples
#'n <- 100
#'family <- "gaussian"
#'y <- stats::rnorm(n=n)
#'y_hat <- stats::rnorm(n=n)
#'calc_metric(y=y,y_hat=y_hat,family=family)
#'
#'family <- "binomial"
#'y <- stats::rbinom(n=n,size=1,prob=0.5)
#'y_hat <- stats::runif(n=n)
#'calc_metric(y=y,y_hat=y_hat,family=family)
#'
calc_metric <- function(y,y_hat,family){
  if(length(y)!=length(y_hat)){
    stop("incompatible lengths")  
  }
  if(family=="gaussian"){
    metric <- mean((y-y_hat)^2)
  } else if(family=="binomial"){
    if(any(y!=0&y!=1)){stop("range")}
    if(any(y_hat<0|y_hat>1)){stop("range")}
    eps <- 1e-06
    metric <- mean(-y*log(pmax(y_hat,eps))-(1-y)*log(1-pmin(y_hat,1-eps)))
  } else {
    stop("family=\"",family,"\" not implemented")
  }
  return(metric)
}

#'@title Create folds for multi-task and transfer learning
#'
#'@export
#'@keywords internal
#'
#'@aliases make_folds_trans
#'
#'@param y
#'multi-task learning:
#'y matrix with \eqn{n} rows (samples) and \eqn{q} columns (outcomes)
#'transfer learning:
#'list of \eqn{q} numeric vectors
#'@param family character \code{"gaussian"} or \code{"binomial"}
#'@param nfolds integer between 2 and \eqn{n}
#'
#'@examples
#'#--- multi-task learning ---
#'family <- "binomial"
#'y <- sim_data_multi(family=family)$y_train
#'fold <- make_folds_multi(y=y,family=family)
#'
#'#--- transfer learning ---
#'family <- "binomial"
#'y <- sim_data_trans(family=family)$y_train
#'fold <- make_folds_trans(y,family=family)
#'
make_folds_multi <- function(y,family,nfolds=10){
  n <- nrow(y)
  q <- ncol(y)
  if(length(family)==1){
    family <- rep(x=family,times=q) 
  }
  if(is.infinite(nfolds)){
    foldid <- sample(seq_len(n))
  } else if(all(family=="gaussian")){
    foldid <- sample(rep(x=seq_len(nfolds),length.out=n))
  } else {
    codes <- apply(y[,family=="binomial",drop=FALSE],1,function(x) paste(x,collapse=""))
    unique <- unique(codes)
    foldid <- rep(x=NA,times=n)
    for(i in unique){
      cond <- codes==i
      if(sum(cond)<nfolds){
        cands <- sample(x=seq_len(nfolds),size=sum(cond))
      } else {
        cands <- seq_len(length.out=nfolds)
      }
      if(sum(cond)==1){ # Should this be sum(cond)<=nfolds?
        foldid[cond] <- cands
      } else {
        foldid[cond] <- sample(x=rep(x=cands,length.out=sum(cond))) # Should this be sample(cands)?
      }
    }
  }
  return(foldid)
}

#'@rdname make_folds_multi
#'@export
#'@keywords internal
#'
make_folds_trans <- function(y,family,nfolds=10){
  q <- length(y)
  if(length(family)==1){
    family <- rep(x=family,times=q) 
  }
  foldid <- list()
  for(i in seq_len(q)){
    if(is.infinite(nfolds)){
      foldid[[i]] <- sample(seq_along(y[[i]]))
    } else {
      if(family[i]=="binomial"){
        foldid[[i]] <- rep(x=NA,times=length(y[[i]]))
        foldid[[i]][y[[i]]==0] <- sample(rep(x=seq_len(nfolds),length.out=sum(y[[i]]==0)))
        foldid[[i]][y[[i]]==1] <- sample(rep(x=seq_len(nfolds),length.out=sum(y[[i]]==1)))
      } else {
        foldid[[i]] <- sample(rep(x=seq_len(nfolds),length.out=length(y[[i]])))
      }
    }
  }
  return(foldid)
}

#'@title 
#'Extract dimensionality.
#'
#'@export
#'@keywords internal
#'
#'@param x list of \eqn{q} matrices,
#'with \eqn{n_k} (samples) rows and \eqn{p} columns (features),
#'for \eqn{k} in \eqn{1,\ldots,q}
#'@param y list of \eqn{q} vectors,
#'with \eqn{n_k} entries,
#'for \eqn{k} in \eqn{1,\ldots,q}
#'
#'@return
#'Returns a list with the slots
#'\eqn{q} (scalar, number of problems),
#'\eqn{n} (vector of length \eqn{q}, number of samples)
#'and \eqn{p} (scalar, number of features)
#'
#'@examples
#'data <- sim_data_trans()
#'get_info(x=data$X_train,y=data$y_train)
#'
get_info <- function(x,y){
  if(length(x)!=length(y)){stop("different q")}
  if(any(sapply(X=x,FUN=base::nrow)!=sapply(X=y,FUN=base::length))){stop("different n")}
  if(any(diff(sapply(X=x,FUN=base::ncol))!=0)){stop("different p")}
  q <- length(x)
  n <- sapply(X=x,FUN=base::nrow)
  p <- ncol(x[[1]])
  list <- list(q=q,n=n,p=p)
  return(list)
}

#'@title Data fusion
#'
#'@export
#'@keywords internal
#'
#'@param x list of \eqn{q} matrices, with \eqn{n_1,\ldots,n_q} rows
#'and with \eqn{p} columns
#'@param y list of \eqn{q} vectors, of length \eqn{n_1,\ldots,n_q},
#'or \code{NULL} (default)
#'@param foldid list of \eqn{q} vectors, of length \eqn{n_1,\ldots,n_q},
#'or \code{NULL} (default)
#'
#'@examples
#'data <- sim_data_trans()
#'sapply(X=data$y_train,FUN=length)
#'sapply(X=data$X_train,FUN=dim)
#'fuse <- fuse_data(x=data$X_train,y=data$y_train)
#'length(fuse$y)
#'dim(fuse$x)
#'table(fuse$index)
#'
fuse_data <- function(x,y=NULL,foldid=NULL){
  list <- list()
  if(!is.null(x)){
    list$x <- do.call(what="rbind",args=x)
  }
  if(!is.null(y)){
    list$y <- do.call(what="c",args=y)
  }
  if(!is.null(foldid)){
    list$foldid <- do.call(what="c",args=foldid)
  }
  q <- length(x)
  n <- sapply(X=x,FUN=base::nrow)
  list$index <- rep(x=seq_len(q),times=n)
  return(list)
}

#'@title
#'Construct internal and external weights
#'
#'@export
#'@keywords internal
#'@param coef matrix with \eqn{p} rows (features) and \eqn{q} columns (problems)
#'@param id integer in \eqn{1,\ldots,q}
#'
#'@return
#'Returns a list with slots
#'\code{lower.limits}, \code{upper.limits},
#'\code{weight.source} (external weights)
#'and \code{weight.target} (internal weights).
#'Each slot is a vector of length \eqn{2*p},
#'with the first \eqn{p} entries for positive effects
#'and the last \eqn{p} entries for negative effects.
#'
#'@seealso
#'Use \code{\link{construct_penfacs}} to obtain penalty factors
#'(i.e., for scaling and inverting weights).
#'
#'@examples
#'p <- 10
#'q <- 3
#'data <- stats::rbinom(p*q,size=1,prob=0.2)*stats::rnorm(p*q)
#'coef <- matrix(data=data,nrow=p,ncol=q)
#'construct_weights(coef=coef,id=1)
#'
construct_weights <- function(coef,id){
  lower.limits <- rep(x=c(0,-Inf),each=nrow(coef))
  upper.limits <- rep(x=c(Inf,0),each=nrow(coef))
  #--- external weights ---
  positive <- apply(X=coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) sum(pmax(0,x)))
  negative <- apply(X=coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) sum(abs(pmin(0,x))))
  weight.ext <- c(positive,negative)
  #--- internal weights ---
  positive <- pmax(0,coef[,id])
  negative <- abs(pmin(0,coef[,id]))
  weight.int <- c(positive,negative)
  list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight.source=weight.ext,weight.target=weight.int)
  return(list)
}

#'@title Construct penalty factors
#'
#'@description
#'Uses internal and external weights
#'as well as internal and external exponents or factors for these weights
#'to construct penalty factors.
#'
#'@details
#'While internal weights are from the problem of interest ("supported" problem),
#'external weights are from the other problems ("supporting" problems).
#'
#'Multiple options exist for scaling prior weights:
#'\itemize{
#'\item \code{"exp"}: \eqn{w_{int}^{v_{int}}+w_{ext}^{v_{ext}}}
#'\item \code{"ari"}: \eqn{v_{int} w_{int} + v_{ext} w_{ext}}
#'\item \code{"geo"}: \eqn{w_{int}^{v_{int}} w_{ext}^{v_{ext}}}
#'\item \code{"rem"}: \eqn{w_{int}^{v_{int}}+w_{ext}^{v_{ext}}-\mathbb{I}(v_{int}=0)-\mathbb{I}(v_{ext}=0))}
#'}
#'The constrained versions \code{"exp.con"}, \code{"ari.con"},
#'\code{"geo.con"}, and \code{"rem.con"} impose \eqn{v_{int}+v_{ext}=1}.
#'The penalty factors are the inverse weights.
#'Suggested choices are \code{"exp"} for predictivity
#'and \code{"ari.con"} for interpretability.
#'
#'@export
#'@keywords internal
#'
#'@param w_int internal weights:
#'numeric vector of length \eqn{p} with non-negative entries
#'@param w_ext external weights:
#'numeric vector of length \eqn{p} with non-negative entries
#'@param v_int exponent or factor for internal weights:
#'non-negative scalar
#'@param v_ext exponent or factor for external weights:
#'non-negative scalar
#'@param type scaling of weights:
#'character \code{"exp"}, \code{"ari"}, \code{"geo"}, or \code{"rem"}
#'(with or without addition of \code{".con"}),
#'default: \code{"exp"}
#'@seealso
#'Use \code{\link{construct_weights}} to obtain \code{w_int} and \code{w_ext}.
#'
#'@examples
#'n <- 10
#'w_int <- stats::runif(n)
#'w_ext <- stats::runif(n)
#'construct_penfacs(w_int,w_ext,v_int=0.5,v_ext=0.5,type="exp")
#'
construct_penfacs <- function(w_int,w_ext,v_int,v_ext,type){
  if(type %in% c("exp","exp.con")){
    pf <- 1/((w_int^v_int)+(w_ext^v_ext))
  } else if(type %in% c("ari","ari.con")){
    pf <- 1/(v_int*w_int + v_ext*w_ext)
  } else if(type %in% c("geo","geo.con")){
    pf <- 1/((w_int^v_int)*(w_ext^v_ext))
  } else if(type %in% c("rem","rem.con")){ 
    pf <- 1/((w_int^v_int)+(w_ext^v_ext)-1*(v_int==0)-1*(v_ext==0))
  }
  if(any(pf<0)){stop("negative pf")}
  return(pf)
}

#----- other methods -----

#'@title
#'Available methods
#'
#'@description
#'Wrapper functions of available methods for related problems.
#'
#'@param x feature matrix (multi-task learning)
#'or list of \eqn{q} feature matrices (transfer learning)
#'@param y response matrix (multi-task learning)
#'or list of \eqn{q} response vectors (transfer learning)
#'@param family character vector with 1 or \eqn{q} entries,
#'possible values are \code{"gaussian"} and sometimes \code{"binomial"} or other
#'@param alpha elastic net mixing parameter:
#'number between 0 and 1
#'@param nfolds number of cross-validation folds: positive integer
#'@param alpha.init elastic net mixing parameter for initial models:
#'number between 0 and 1
#'@param lambda sequence of regularisation parameters
#'@param object output from multi-task learning or transfer learning method
#'@param newx feature matrix (MTL) or list of feature matrices (TL)
#'of testing samples
#'@param ... (not applicable)
#'
#'@return
#'The wrapper functions \code{wrap_empty}, \code{wrap_separate},
#'\code{wrap_common}, \code{wrap_mgaussian}, \code{wrap_spls},
#'\code{wrap_glmtrans}, and \code{wrap_xrnet} return fitted models,
#'and the generic functions \code{coef} and \code{predict}
#'return coefficients or predicted values in a standardised format.
#'
#'@examples
#'#--- multi-task learning ---
#'n_train <- 100
#'n_test <- 10
#'p <- 50
#'q <- 3
#'family <- "gaussian"
#'x <- matrix(data=rnorm(n=n_train*p),nrow=n_train,ncol=p)
#'newx <- matrix(data=rnorm(n=n_test*p),nrow=n_test,ncol=p)
#'y <- matrix(data=rnorm(n_train*q),nrow=n_train,ncol=q)
#'object <- wrap_empty(x=x,y=y,family=family)
#'model <- "empty" # try "empty", "separate", "mgaussian" or "spls"
#'if(model=="empty"){
#'   object <- wrap_empty(x=x,y=y,family=family)
#'} else if(model=="separate"){
#'   object <- wrap_separate(x=x,y=y,family=family)
#'} else if(model=="mgaussian"){
#'   object <- wrap_mgaussian(x=x,y=y,family=family)
#'} else if(model=="spls"){
#'   object <- wrap_spls(x=x,y=y,family=family)
#'}
#'coef(object)
#'predict(object,newx=newx)
#'
#'#--- transfer learning ---
#'n_train <- c(100,50)
#'n_test <- c(10,10)
#'p <- 50
#'x <- lapply(X=n_train,function(n) matrix(data=stats::rnorm(n*p),nrow=n,ncol=p))
#'newx <- lapply(X=n_test,function(n) matrix(data=stats::rnorm(n*p),nrow=n,ncol=p))
#'y <- lapply(X=n_train,function(n) stats::rnorm(n))
#'family <- "gaussian"
#'model <- "empty" # try "empty", "separate", "common", "glmtrans", or "xrnet"
#'if(model=="empty"){
#'  object <- wrap_empty(x=x,y=y,family=family)
#'} else if(model=="separate"){
#'  object <- wrap_separate(x=x,y=y,family=family)
#'} else if(model=="common"){
#'  object <- wrap_common(x=x,y=y,family=family)
#'} else if(model=="glmtrans"){
#'  object <- wrap_glmtrans(x=x,y=y,family=family)
#'} else if(model=="xrnet"){
#'  object <- wrap_xrnet(x=x,y=y,family=family)
#'}
#'coef(object)
#'predict(object,newx=newx)
#'
#'@seealso
#'See original functions
#'\link[glmnet]{cv.glmnet} (with argument \code{family="mgaussian"})
#'\link[spls]{spls}, \link[glmtrans]{glmtrans}, and \link[xrnet]{xrnet}.
#'
#'@references
#'Noah Simon, Jerome H. Friedman, and Trevor Hastie (2013). 
#"A blockwise descent algorithm for group-penalized multiresponse and multinomial regression".
#'\emph{arXiv} (Preprint).
#'\doi{10.48550/arXiv.1311.6529}.
#'(\link[glmnet]{cv.glmnet})
#'
#'Hyonho Chun and Sündüz Keleş (2010).
#'"Sparse Partial Least Squares Regression for Simultaneous Dimension Reduction and Variable Selection".
#'\emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}
#'72(1);3–25.
#'\doi{10.1111/j.1467-9868.2009.00723.x}.
#'(\link[spls]{spls})
#'
#'Ye Tian and Yang Feng (2022).
#'"Transfer learning under high-dimensional generalized linear models".
#'\emph{Journal of the American Statistical Association}
#'118(544):2684-2697.
#'\doi{10.1080/01621459.2022.2071278}.
#'(\link[glmtrans]{glmtrans})
#'
#'Garrett M. Weaver and Juan Pablo Lewinger (2019).
#'"xrnet: Hierarchical Regularized Regression to Incorporate External Data".
#'\emph{Journal of Open Source Software}
#'4(44):1761.
#'\doi{10.21105/joss.01761}.
#'(\link[xrnet]{xrnet})
#'
#'@name methods
#'
NULL

#'@describeIn methods intercept-only model (MTL and TL)
#'@export
#'@keywords internal
wrap_empty <- function(x,y,family,alpha=1){
  object <- wrap_separate(x=x,y=y,family=family,alpha=alpha,lambda=c(99e99,99e98))
  return(object)
}

#'@describeIn methods separate model for each problem (MTL and TL)
#'@export
#'@keywords internal
#'
wrap_separate <- function(x,y,family,alpha=1,lambda=NULL){
  if(is.matrix(x) & is.matrix(y)){
    q <- ncol(y)
    x <- replicate(n=q,expr=x,simplify=FALSE)
    y <- apply(y,2,function(x) x,simplify=FALSE)
  }
  if(length(family)==1){
    family <- rep(x=family,times=length(y))
  }
  object <- list()
  object$info <- get_info(x=x,y=y)
  object$cv.glmnet <- list()
  for(i in seq_len(object$info$q)){
    object$cv.glmnet[[i]] <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],family=family[i],alpha=alpha,lambda=lambda)
  }
  class(object) <- "wrap_separate"
  return(object)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
predict.wrap_separate <- function(object,newx,...){
  if(is.matrix(newx)){
    newx <- replicate(n=object$info$q,expr=newx,simplify=FALSE) 
  }
  q <- length(newx)
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- stats::predict(object$cv.glmnet[[i]],newx=newx[[i]],s="lambda.min",type="response")
  }
  return(y_hat)
}

#'@rdname methods
#'@export
#'@keywords internal
coef.wrap_separate <- function(object,...){
  p <- object$info$p
  q <- object$info$q
  alpha <- rep(x=NA,times=q)
  beta <- matrix(data=NA,nrow=p,ncol=q)
  for(i in seq_len(q)){
    coef <- stats::coef(object$cv.glmnet[[i]],s="lambda.min")
    alpha[i] <- coef[1]
    beta[,i] <- coef[-1]
  }
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

#'@describeIn methods common model for all problems (TL)
#'@export
#'@keywords internal
#'
wrap_common <- function(x,y,family,alpha=1){
  family <- unique(family)
  if(length(family)>1){stop("requires unique family")}
  object <- list()
  object$info <- get_info(x=x,y=y)
  fuse <- fuse_data(x=x,y=y,foldid=NULL)
  object$cv.glmnet <- glmnet::cv.glmnet(x=fuse$x,y=fuse$y,family=family,alpha=alpha)
  class(object) <- "wrap_common"
  return(object)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
predict.wrap_common <- function(object,newx,...){
  fuse <- fuse_data(x=newx,y=NULL,foldid=NULL)
  temp <- stats::predict(object=object$cv.glmnet,newx=fuse$x,s=object$cv.glmnet$lambda.min,type="response")
  y_hat <- tapply(X=temp,INDEX=fuse$index,FUN=function(x) x)
  return(y_hat)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
coef.wrap_common <- function(object,...){
  coef <- stats::coef(object=object$cv.glmnet,s="lambda.min")
  alpha <- rep(x=coef[1],times=object$info$q)
  beta <- matrix(data=coef[-1],nrow=object$info$p,ncol=object$info$q)
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

#'@describeIn methods multivariate Gaussian regression (MTL)
#'@export
#'@keywords internal
#'
wrap_mgaussian <- function(x,y,family="gaussian",alpha=1){
  object <- list()
  family <- unique(family)
  if(length(family)>1){stop("requires unique family")}
  if(family!="gaussian"){stop("requires gaussian")}
  object$cv.glmnet <- glmnet::cv.glmnet(x=x,y=y,family="mgaussian",alpha=alpha)
  class(object) <- "wrap_mgaussian"
  return(object)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
predict.wrap_mgaussian <- function(object,newx,...){
  y_hat <- stats::predict(object$cv.glmnet,newx=newx,s="lambda.min")
  apply(y_hat,2,function(x) x,simplify=FALSE)
}

#'@rdname methods
#'@export
#'@keywords internal
coef.wrap_mgaussian <- function(object,...){
  coef <- stats::coef(object$cv.glmnet,s="lambda.min")
  alpha <- sapply(coef,function(x) x[1])
  beta <- sapply(coef,function(x) x[-1])
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

#'@describeIn methods sparse partial least squares (MTL)
#'@export
#'@keywords internal
#'
wrap_spls <- function(x,y,family="gaussian",alpha=1,nfolds=10){
  if(any(family!="gaussian")){stop("SPLS requires Gaussian family.")}
  if(alpha!=1){stop("SPLS requires alpha=1")}
  invisible(utils::capture.output(cv.spls <- spls::cv.spls(x=x,y=y,fold=10,K=seq_len(min(ncol(x),10)),
                                                           eta=seq(from=0.0,to=0.9,by=0.1),plot.it=FALSE)))
  object <- spls::spls(x=x,y=y,K=cv.spls$K.opt,eta=cv.spls$eta.opt)
  class(object) <- "wrap_spls"
  return(object)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
predict.wrap_spls <- function(object,newx,...){
  temp <- spls::predict.spls(object=object,newx=newx,type="fit")
  y_hat <- apply(X=temp,MARGIN=2,FUN=function(x) x,simplify=FALSE)
  return(y_hat)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
coef.wrap_spls <- function(object,...){
  coef <- spls::coef.spls(object=object)
  list <- list(alpha=NA,beta=coef)
  return(list)
}

#'@describeIn methods transfer generalised linear model (TL)
#'@export
#'@keywords internal
#'
wrap_glmtrans <- function(x,y,family="gaussian",alpha=1){
  family <- unique(family)
  if(length(family)>1){stop("glmtrans-wrapper requires unique family")}
  q <- length(x)
  source <- list()
  for(i in seq_len(q)){
    source[[i]] <- list(y=y[[i]],x=x[[i]])
  }
  object <- list()
  for(i in seq_len(q)){
    invisible(utils::capture.output(object[[i]] <- glmtrans::glmtrans(target=list(y=y[[i]],x=x[[i]]),source=source[-i],family=family,alpha=alpha)))
  }
  class(object) <- "wrap_glmtrans"
  return(object)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
predict.wrap_glmtrans <- function(object,newx,...){
  q <- length(newx)
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- stats::predict(object=object[[i]],newx=newx[[i]],type="response",s="lambda.min")
  }
  return(y_hat)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
coef.wrap_glmtrans <- function(object,...){
  coef <- sapply(object,function(x) x$beta)
  alpha <- coef[1,]
  beta <- coef[-1,]
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

#'@describeIn methods hierarchical regression (TL) 
#'@export
#'@keywords internal
#'
wrap_xrnet <- function(x,y,alpha.init=0.95,alpha=1,nfolds=10,family="gaussian"){
  family <- unique(family)
  if(length(family)!=1){stop("xrnet-wrapper requires single family")}
  q <- length(x)
  p <- ncol(x[[1]])
  #--- stage 1 ---
  prior <- matrix(data=NA,nrow=p,ncol=q)
  for(i in seq_len(q)){
    object <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],alpha=alpha.init,family=family)
    prior[,i] <- stats::coef(object=object,s="lambda.min")[-1]
  }
  #--- stage 2 ---
  object <- list()
  for(i in seq_len(q)){
    #if(family=="binomial"){
    #  foldid <- rep(x=NA,times=length(y[[i]]))
    #  foldid[y[[i]]==0] <- sample(rep(x=seq_len(nfolds),length.out=sum(y[[i]]==0)))
    #  foldid[y[[i]]==1] <- sample(rep(x=seq_len(nfolds),length.out=sum(y[[i]]==1)))
    #} else {
    #  foldid <- NULL
    #}
    # use coefficients from other problems as prior information
    cond <- apply(prior,2,stats::sd)>0 & seq_len(q)!=i
    for(j in 1:2){
      temp <- tryCatch(expr=xrnet::tune_xrnet(x=x[[i]],y=y[[i]],
                                  external=prior[,cond,drop=FALSE],
                                  penalty_main=xrnet::define_penalty(penalty_type=alpha),
                                  family=family,
                                  #foldid=foldid, # was NULL
                                  standardize=c(j==1,j==1),
                                  nfolds=nfolds,
                                  loss="deviance"),error=function(x) NULL)
      #if(is.null(temp)){
      #  cat("Run",j,"returns NULL.\n")
      #}
      if(!is.null(temp)){break}
    }
    object[[i]] <- temp
  }
  class(object) <- "wrap_xrnet"
  return(object)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
predict.wrap_xrnet <- function(object,newx,...){
  y_hat <- list()
  for(i in seq_along(object)){
    y_hat[[i]] <- stats::predict(object[[i]],newdata=newx[[i]])
  }
  return(y_hat)
}

#'@rdname methods
#'@export
#'@keywords internal
#'
coef.wrap_xrnet <- function(object,...){
  alpha <- beta <- numeric()
  for(i in seq_along(object)){
    coef <- stats::coef(object[[i]])
    alpha[i] <- coef$beta0
    beta <- cbind(beta,coef$betas)
  }
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

#----- simulation and application -----

#'@title Data simulation for related problems
#'
#'@description
#'Simulates data for multi-task learning and transfer learning.
#'
#'@export
#'@keywords internal
#'
#'@aliases sim_data_trans
#'
#'@param prob.common probability of common effect (number between 0 and 1)
#'@param prob.separate probability of separate effect (number between 0 and 1)
#'@param q number of datasets: integer
#'@param n0 number of training samples: integer vector of length \eqn{q}
#'@param n1 number of testing samples for all datasets: integer
#'@param p number of features: integer
#'@param family character \code{"gaussian"} or \code{"binomial"}
#'@param rho correlation (for decreasing structure)
#'
#'@return
#'\itemize{
#'\item Multi-task learning:
#'Returns a list with slots
#'\code{y_train} (\eqn{n_0 \times q} matrix),
#'\code{X_train}(\eqn{n_0 \times p} matrix),
#'\code{y_test} (\eqn{n_1 \times q} matrix),
#'\code{X_test} (\eqn{n_1 \times p} matrix),
#'and \code{beta} (\eqn{p \times q} matrix).
#'\item Transfer learning:
#'Returns a list with slots
#'\code{y_train} (\eqn{q} vectors)
#'and \code{X_train} (\eqn{q} matrices with \eqn{p} columns) for training data,
#'and \code{y_test} (\eqn{vectors}) and \code{X_test}
#'(\eqn{q} matrices with \eqn{p} columns)
#'for testing data,
#'and \code{beta} for effects (\eqn{p \times q} matrix).
#'}
#'
#'@examples
#'#--- multi-task learning ---
#'data <- sim_data_multi()
#'sapply(X=data,FUN=dim)
#'
#'#--- transfer learning ---
#'data <- sim_data_trans()
#'sapply(X=data$y_train,FUN=length)
#'sapply(X=data$X_train,FUN=dim)
#'sapply(X=data$y_test,FUN=length)
#'sapply(X=data$X_test,FUN=dim)
#'dim(data$beta)
#'
sim_data_multi <- function(prob.common=0.05,prob.separate=0.05,q=3,n0=100,n1=10000,p=200,rho=0.5,family="gaussian"){
  n <- n0 + n1
  theta <- stats::rnorm(n=p)*stats::rbinom(n=p,size=1,prob=prob.common)
  if(rho==0){
    X <- matrix(data=stats::rnorm(n*p),nrow=n,ncol=p)
  } else {
    mean <- rep(x=0,times=p)
    sigma <- matrix(data=NA,nrow=p,ncol=p)
    sigma <- rho^abs(row(sigma)-col(sigma))
    X <- mvtnorm::rmvnorm(n=n,mean=mean,sigma=sigma)
  }
  y <- matrix(data=NA,nrow=n,ncol=q)
  beta <- matrix(data=NA,nrow=p,ncol=q)
  for(i in seq_len(q)){
    beta[,i] <- theta + stats::rnorm(n=p)*stats::rbinom(n=p,size=1,prob=prob.separate)
    y[,i] <- X %*% beta[,i] + stats::rnorm(n=n)
    if(family=="binomial"){
      y[,i] <- 1*(y[,i]>=0)
    }
  }
  foldid <- rep(x=c(0,1),times=c(n0,n1))
  cond <- foldid==0
  y_train <- y[cond,]
  y_test <- y[!cond,]
  X_train <- X[cond,]
  X_test <- X[!cond,]
  list <- list(y_train=y_train,X_train=X_train,y_test=y_test,X_test=X_test,beta=beta)
  return(list)
}

#'@rdname sim_data_multi
#'@export
#'@keywords internal
#'
sim_data_trans <- function(prob.common=0.05,prob.separate=0.05,q=3,n0=c(50,100,200),n1=10000,p=200,rho=0.5,family="gaussian"){
  if(length(n0)==1){
    n0 <- rep(x=n0,times=q)
  } else {
    if(length(n0)!=q){
      stop("Invalid.")
    }
  }
  n <- n0 + n1
  theta <- stats::rnorm(n=p)*stats::rbinom(n=p,size=1,prob=prob.common)
  X <- beta <- y <- foldid <- list()
  beta <- matrix(data=NA,nrow=p,ncol=q)
  for(i in seq_len(q)){
    if(rho==0){
      X[[i]] <- matrix(data=stats::rnorm(n[i]*p),nrow=n[i],ncol=p)
    } else {
      mean <- rep(x=0,times=p)
      sigma <- matrix(data=NA,nrow=p,ncol=p)
      sigma <- rho^abs(row(sigma)-col(sigma))
      X[[i]] <- mvtnorm::rmvnorm(n=n[i],mean=mean,sigma=sigma)
    }
    beta[,i] <- theta + stats::rnorm(n=p)*stats::rbinom(n=p,size=1,prob=prob.separate)
    y[[i]] <- X[[i]] %*% beta[,i] + stats::rnorm(n=n[[i]])
    if(family=="binomial"){
      y[[i]] <- 1*(y[[i]]>=0)
    }
    foldid[[i]] <- rep(x=c(0,1),times=c(n0[i],n1))
  }
  y_train <- X_train <- y_test <- X_test <- list()
  for(i in seq_len(q)){
    cond <- foldid[[i]]==0
    y_train[[i]] <- y[[i]][cond]
    y_test[[i]] <- y[[i]][!cond]
    X_train[[i]] <- X[[i]][cond,]
    X_test[[i]] <- X[[i]][!cond,]
  }
  list <- list(y_train=y_train,X_train=X_train,y_test=y_test,X_test=X_test,beta=beta)
  return(list)
}

#'@title Train and test model
#'
#'@description
#'Trains and tests prediction models
#'
#'@export
#'@keywords internal
#'
#'@param y_train target of training samples:
#'\eqn{n \times q} matrix (multi-task learning)
#'or list of \eqn{q} vectors of length \eqn{n_1,\ldots,n_q} (transfer learning)
#'@param X_train features of training samples:
#'\eqn{n \times p} matrix (multi-task learning)
#'or list of \eqn{q} matrices of dimensions
#'\eqn{n_1 \times p,\ldots,n_q \times p} (transfer learning)
#'@param y_test target of testing samples:
#'\eqn{m \times p} matrix (multi-task learning)
#' or list of \eqn{q} vectors of length \eqn{m_1,\ldots,m_q} (transfer learning)
#'@param X_test features of testing samples:
#' \eqn{m \times p} matrix (multi-task learning) or
#'list of \eqn{q} matrices of dimensions
#'\eqn{m_1 \times p,\ldots,m_q \times p} (transfer learning)
#'@inheritParams sparselink
#'
#'@examples
#'#--- multi-task learning ---
#'\donttest{
#'family <- "gaussian"
#'data <- sim_data_multi(family=family)
#'result <- traintest(data$y_train,data$X_train,family=family)
#'}
#'
#'#--- transfer learning ---
#'\donttest{
#'family <- "gaussian"
#'data <- sim_data_trans(family=family)
#'result <- traintest(data$y_train,data$X_train,family=family)
#'}
#'
traintest <- function(y_train,X_train,y_test=NULL,X_test=NULL,family="gaussian",alpha=1,method=c("wrap_empty","wrap_separate","sparselink"),alpha.init=0.95,type="exp",cands=NULL){
  if(is.list(y_train)){
    q <- length(y_train)
  } else {
    q <- ncol(y_train)
  }
  if(length(family)==1){
    family <- rep(x=family,times=q)
  }
  if(is.matrix(y_test)){
    y_test <- apply(y_test,2,function(x) x,simplify=FALSE)
  }
  time <- rep(x=0,times=length(method))
  names(time) <- method
  deviance <- auc <- matrix(data=NA,nrow=length(method),ncol=q,dimnames=list(method,NULL))
  coef <- y_hat <- list()
  for(i in seq_along(method)){
    cat("method",method[i],"\n")
    func <- eval(parse(text=paste0(method[i])))
    start <- Sys.time()
    hyperpar <- NULL
    if(method[i]=="sparselink"){
      object <- func(x=X_train,y=y_train,family=family,alpha.init=alpha.init,alpha=alpha,type=type,cands=cands)
      hyperpar <- object$weight.min
    } else if(method[i]=="devel"){
      object <- func(x=X_train,y=y_train,family=family,alpha.init=alpha.init,alpha=alpha)
    } else {
      object <- func(x=X_train,y=y_train,family=family,alpha=alpha)
    }
    end <- Sys.time()
    time[i] <- difftime(time1=end,time2=start,units="secs")
    if(!is.null(X_test)){
      y_hat[[i]] <- stats::predict(object,newx=X_test)
    }
    if(is.null(y_test)){
      deviance[i,] <- auc[i,] <- NA
    } else {
      for(j in seq_len(q)){
        deviance[i,j] <- calc_metric(y=y_test[[j]],y_hat=y_hat[[i]][[j]],family=family[j])
        if(family[j]=="binomial"){
          auc[i,j] <- pROC::auc(response=y_test[[j]],predictor=as.vector(y_hat[[i]][[j]]),direction="<",levels=c(0,1))
        }
      }
    }
    coef[[i]] <- stats::coef(object)$beta
  }
  names(coef) <- method
  if(!is.null(X_test)){
    names(y_hat) <- method
  }
  list <- list(time=time,deviance=deviance,auc=auc,coef=coef,y_hat=y_hat,hyperpar=hyperpar)
  return(list)
}

#'@title Model comparison
#'
#'@description
#'Compares predictive methods for multi-task learning (\code{cv_multiple}) or
#'transfer learning (\code{cv_transfer}) by \eqn{k}-fold cross-validation.
#'
#'@export
#'@keywords internal
#'
#'@aliases cv_transfer
#'
#'@inheritParams sparselink
#'
#'@examples
#'#--- multi-task learning ---
#'\donttest{
#'family <- "gaussian"
#'data <- sim_data_multi(family=family)
#'metric <- cv_multiple(y=data$y_train,X=data$X_train,family=family)
#'metric$deviance}
#'
#'#--- transfer learning ---
#'\donttest{
#'family <- "gaussian"
#'data <- sim_data_trans(family=family)
#'metric <- cv_transfer(y=data$y_train,X=data$X_train,family=family)
#'metric$deviance}
#'
cv_multiple <- function(y,X,family,alpha=1,nfolds=10,method=c("wrap_separate","wrap_mgaussian","sparselink","wrap_spls"),alpha.init=0.95,type="exp",cands=NULL){
  mode <- "multiple"
  foldid <- make_folds_multi(y=y,family=family,nfolds=nfolds)
  n <- nrow(y)
  q <- ncol(y)
  
  y_hat <- list()
  for(k in method){
    y_hat[[k]] <- matrix(data=NA,nrow=nrow(y),ncol=q)
  }
  
  for(i in seq_len(nfolds)){
    cat("fold",i,"\n")
    y_train <- y[foldid!=i,]
    y_test <- y[foldid==i,]
    X_train <- X[foldid!=i,,drop=FALSE]
    X_test <- X[foldid==i,,drop=FALSE]
    test <- traintest(y_train=y_train,X_train=X_train,X_test=X_test,y_test=y_test,family=family,method=method,alpha=alpha,alpha.init=alpha.init,type=type,cands=cands)
    for(k in method){
      y_hat[[k]][foldid==i,] <- as.matrix(as.data.frame(test$y_hat[[k]]))
    }
  }
  
  deviance <- auc <- matrix(data=NA,nrow=ncol(y),ncol=length(method),dimnames=list(names(y),method))
  for(j in seq_len(ncol(y))){
    for(k in method){
      deviance[j,k] <- calc_metric(y=y[,j],y_hat=y_hat[[k]][,j],family=family)
      if(family=="binomial"){
        auc[j,k] <- pROC::auc(response=y[,j],predictor=as.vector(y_hat[[k]][,j]),direction="<",levels=c(0,1))
      }
    }
  }
  # refit model on all folds
  #cat("refit on all folds","\n")
  #refit <- traintest(y_train=y,X_train=X,family=family,method=method,alpha.init=alpha.init,type=type,alpha=alpha)
  list <- list(deviance=deviance,auc=auc)
  return(list)
}

#'@rdname cv_multiple
#'@export
#'@keywords internal
#'
cv_transfer <- function(y,X,family,alpha=1,nfolds=10,method=c("wrap_separate","wrap_glmtrans","sparselink","wrap_xrnet"),alpha.init=0.95,type="exp",cands=NULL){
  mode <- "transfer"
  foldid <- make_folds_trans(y=y,family=family,nfolds=nfolds)
  n <- length(y[[1]])
  q <- length(y)
  
  y_hat <- list()
  for(j in seq_len(q)){
    y_hat[[j]] <- matrix(data=NA,nrow=length(y[[j]]),ncol=length(method),dimnames=list(NULL,method))
  }
  for(i in seq_len(nfolds)){
    cat("fold",i,"\n")
    y_train <- X_train <- X_test <- list()
    for(j in seq_len(q)){
      y_train[[j]] <- y[[j]][foldid[[j]]!=i]
      X_train[[j]] <- X[[j]][foldid[[j]]!=i,,drop=FALSE]
      X_test[[j]] <- X[[j]][foldid[[j]]==i,,drop=FALSE]
    }
    test <- traintest(y_train=y_train,X_train=X_train,X_test=X_test,family=family,method=method,alpha=alpha,alpha.init=alpha.init,type=type,cands=cands)
    for(j in seq_len(q)){
      for(k in method){
        y_hat[[j]][foldid[[j]]==i,k] <- test$y_hat[[k]][[j]]
      }
    }
  }
  deviance <- auc <- matrix(data=NA,nrow=length(y),ncol=length(method),dimnames=list(names(y),method))
  for(j in seq_along(y)){
    for(k in method){
      deviance[j,k] <- calc_metric(y=y[[j]],y_hat=y_hat[[j]][,k],family=family)
      if(family=="binomial"){
        auc[j,k] <- pROC::auc(response=y[[j]],predictor=as.vector(y_hat[[j]][,k]),direction="<",levels=c(0,1))
      }
    }
  }
  # refit model on all folds
  cat("refit on all folds","\n")
  refit <- traintest(y_train=y,X_train=X,family=family,method=method,alpha.init=alpha.init,type=type,alpha=alpha,cands=cands)
  list <- list(deviance=deviance,auc=auc,refit=refit)
  return(list)
}

#'@title Metrics for sign detection
#'
#'@description
#'Calculates sensitivity, specificity and precision for ternary data
#'(with -1 for negative effect, 0 for no effect, 1 for positive effect).
#'
#'@export
#'@keywords internal
#'
#'@aliases count_matrix
#'
#'@param truth (i) vector of length \eqn{p} or
#'(ii) \eqn{n \times p} matrix with entries in -1, 0, 1
#'@param estim (i) vector of length \eqn{p} or
#'(ii) \eqn{n \times p} matrix with entries -1, 0, 1
#'
#'@examples
#'truth <- sample(x=c(-1,0,1),size=20,replace=TRUE)
#'estim <- sample(x=c(-1,0,1),size=20,replace=TRUE)
#'table(truth,estim)
#'count_vector(truth,estim)
#'
count_vector <- function(truth,estim){
  if(!is.vector(truth)){stop()}
  if(length(truth)!=length(estim)){stop()}
  if(!all(truth %in% c(-1,0,1))){stop()}
  if(!all(estim %in% c(-1,0,1))){stop()}
  TN <- sum(truth==0 & estim==0)
  TP <- sum(truth!=0 & estim==truth)
  #FN <- sum(truth!=0 & estim==0) # original
  FN <- sum(truth!=0 & estim!=truth)
  FP <- sum(estim!=0 & estim!=truth)
  DD <- sum(truth!=0 & truth==-estim)
  if(TN+TP+FN+FP-DD!=length(truth)){stop()}
  sensitivity <- (sum(truth==1 & estim==1)+sum(truth==-1 & estim==-1))/sum(truth==1|truth==-1)
  specificity <- sum(truth==0 & estim==0)/sum(truth==0)
  #return(c(TN=TN,TP=TP,FN=FN,FP=FP)/length(truth))
  precision <- (sum(truth==1 & estim==1)+sum(truth==-1 & estim==-1))/(sum(estim==1 | estim==-1))
  return(c(sensitivity=sensitivity,specificity=specificity,precision=precision))
}

#'@rdname count_vector
#'@export
#'@keywords internal
#'
count_matrix <- function(truth,estim){
  if(!is.matrix(truth)){stop()}
  if(any(dim(truth)!=dim(estim))){stop()}
  rate <- numeric()
  for(i in 1:ncol(truth)){
    rate <- cbind(rate,count_vector(truth=truth[,i],estim=estim[,i]))
  }
  return(rate)
}

#----- graphics -----

#'@title Pairwise differences
#'
#'@description
#'Visualises differences within sets of three values in different settings.
#'
#'@export
#'@keywords internal
#'
#'@param x setting: character vector
#'@param y0 values on the left: numeric vector
#'@param y1 values in the centre: numeric vector
#'@param y2 values on the right: numeric vector
#'@param dist horizontal distance between points
#'@param main title
#'@param increase change to arrow NULL, up, down
#'@param cex.axis numeric
#'@param cex.main numeric
#'
#'@examples
#'m <- 3 # number of settings
#'n <- 5 # number of repetitions
#'x <- rep(LETTERS[1:m],each=n)
#'y0 <- stats::rnorm(n*m,mean=0)
#'y1 <- stats::rnorm(n*m,mean=ifelse(x=="A",2,-2))
#'y2 <- stats::rnorm(n*m,mean=ifelse(x=="A",4,-4))
#'plot_change(x,y0,y1,y2)
#'
plot_change <- function(x,y0,y1,y2,dist=0.15,main="",cex.axis=0.5,cex.main=1,increase=TRUE){
  unique <- unique(x)
  graphics::plot.new()
  xlim <- c(1-0.2,length(unique)+0.2)
  ylim <- range(c(y0,y1,y2),na.rm=TRUE)
  graphics::plot.window(xlim=xlim,ylim=ylim)
  graphics::mtext(text=unique,side=1,at=seq_along(unique),line=1,cex=cex.axis)  
  graphics::axis(side=2,cex.axis=cex.axis)
  for(i in seq_along(unique)){
    cond <- x==unique[i]
    graphics::segments(x0=i,y0=y1[cond],x1=i+dist,y1=y2[cond],col="grey")
    graphics::segments(x0=i-dist,y0=y0[cond],x1=i,y1=y1[cond],col="grey")
    graphics::points(x=rep(i-dist,times=sum(cond)),y=y0[cond],col="red",pch=16,cex=0.8)
    graphics::points(x=rep(i,times=sum(cond)),y=y1[cond],col="blue",pch=16,cex=0.8)
    graphics::points(x=rep(i+dist,times=sum(cond)),y=y2[cond],col="red",pch=16,cex=0.8)
  }
  graphics::title(main=main,line=0,cex.main=cex.main)
  graphics::par(xpd=TRUE)
  usr <- graphics::par("usr")
  margin <- 0.1*diff(usr[3:4])
  if(increase){
    inferior <- usr[3]
    superior <- usr[4]
    margin.inferior <- +margin
    margin.superior <- -margin
  } else {
    inferior <- usr[4]
    superior <- usr[3]
    margin.inferior <- -margin
    margin.superior <- +margin
  }
  pos <- xlim[1]-0.13*diff(xlim)
  pos <- usr[2] # trial
  graphics::arrows(x0=pos,y0=inferior+2*margin.inferior,
                   y1=superior+2*margin.superior,lwd=2,length=0.08,col="grey")
  graphics::text(x=pos,y=inferior+1*margin.inferior,labels="-",col="red",font=2,cex=1.2)
  graphics::text(x=pos,y=superior+1*margin.superior,labels="+",col="blue",font=2,cex=1.2)
  graphics::par(xpd=FALSE)
}

#'@title Visualise metric that depends on two parameters
#'@export
#'@keywords internal
#'
#'@description 
#'Displays values in \code{y} in a grey scale (white=lowest, black=highest),
#'for different combinations of the two variables in \code{x}.
#'The lowest value is indicated by a red cross,
#'and the lowest value on the diagonal is indicated by a red circle.
#'
#'@param x list with slots source and target
#'@param y numeric vector
#'
#'@examples
#'values <- seq(from=0,to=1,by=0.2)
#'x <- expand.grid(source=values,target=values)
#'y <- stats::rexp(n=length(values)*length(values))
#'plot_weight(x=x,y=y)
#'
plot_weight <- function(x,y){
  if(stats::cor(x$source,x$target)==-1){
    graphics::plot(x=x$source,y=y,type="o",xlab="weight source = 1 - weight target")
  } else {
    col <- grDevices::grey(level=1-(y-min(y))/(max(y)-min(y)),alpha=1)
    graphics::plot.new()
    graphics::plot.window(xlim=range(x$source),ylim=range(x$target))
    graphics::box()
    graphics::axis(side=1)
    graphics::axis(side=2)
    graphics::title(xlab="weight source",ylab="weight target")
    graphics::abline(a=1,b=-1)
    for(i in seq_len(nrow(x))){
      graphics::points(x=x$source[i],y=x$target[i],col="black",bg=col[i],pch=21,cex=3)
    }
    id <- which.min(x=y)
    graphics::points(x=x$source[id],y=x$target[id],pch=4,col="red",lwd=2,cex=2)
    diagonal <- x$source+x$target==1
    id <- which(y==min(y[diagonal]) & diagonal)
    graphics::points(x=x$source[id],y=x$target[id],pch=1,col="red",lwd=2,cex=1)
  }
}
