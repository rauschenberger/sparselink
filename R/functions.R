
if(FALSE){
  setwd("C:/Users/arauschenberger/Desktop/sparselink/package")
  roxygen2::roxygenise()
  rcmdcheck::rcmdcheck()
}

#' @title logit function
#' @export
#' @keywords internal
#' 
#' @param x numeric vector with values in unit interval
#' 
#' @examples
#' x <- seq(from=0,to=1,length.out=100)
#' y <- logit(x=x)
#' graphics::plot(x=x,y=y,type="l")
#' graphics::abline(v=0.5,lty=2)
#' graphics::abline(h=0,lty=2)
#' 
logit <- function(x){
  if(any(x<0|x>1)){stop("x must be in unit interval")}
  log(x/(1-x))
}

#' @title Sigmoid function
#' @export
#' @keywords internal
#' 
#' @param x numeric vector
#'
#' @examples
#' x <- seq(from=-3,to=3,length.out=100)
#' y <- sigmoid(x)
#' graphics::plot(x=x,y=y,type="l")
#' graphics::abline(v=0,lty=2)
#' graphics::abline(h=0.5,lty=2)
#' 
sigmoid <- function(x){
  1/(1+exp(-x))
}

#' @title Link function
#' @export
#' @keywords internal
#' 
#' @description
#' Applies the link function.
#'
#' @param mu numeric vector (with values in unit interval if family="binomial")
#' @param family character "gaussian" or "binomial"
#'
#' @examples
#' family <- "binomial"
#' from <- ifelse(family=="binomial",0,-3)
#' to <- ifelse(family=="binomial",1,3)
#' mu <- seq(from=from,to=to,length.out=100)
#' eta <- link_function(mu=mu,family=family)
#' graphics::plot(x=mu,y=eta,type="l",main=family)
#' v <- ifelse(family=="binomial",0.5,0)
#' graphics::abline(v=v,lty=2)
#' graphics::abline(h=0,lty=2)
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

#' @title Mean function
#' @export
#' @keywords internal
#' 
#' @description
#' Applies the mean function (inverse link function).
#' 
#' @param eta numeric vector
#' @param family character "gaussian" or "binomial"
#' 
#' @examples
#' family <- "binomial"
#' eta <- seq(from=-3,to=3,length.out=100)
#' mu <- mean_function(eta=eta,family=family)
#' graphics::plot(x=eta,y=mu,type="l",main=family)
#' graphics::abline(v=0,lty=2)
#' h <- ifelse(family=="binomial",0.5,0)
#' graphics::abline(h=h,lty=2)
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

#' @title Data simulation for transfer learning
#' 
#' @description
#' Simulates data for transfer learning.
#' 
#' @export
#' 
#' @param prob.common probability of common effect
#' @param prob.separate probability of separate effect
#' @param q number of datasets: integer
#' @param n0 number of training samples: integer vector of length q
#' @param n1 number of testing samples for all datasets: integer
#' @param p number of features: integer
#' @param family character "gaussian" or "binomial"
#' @param rho correlation (for decreasing structure)
#' 
#' @returns
#' Returns a list with slots y_train and X_train for training data,
#' y_test and X_test for testing data,
#' and beta for effects.
#' The training data contains vectors of different lengths (y_train)
#' and matrices with different number of rows (X_train).
#' 
#' @examples
#' data <- sim.data.transfer()
#' sapply(X=data$y_train,FUN=length)
#' sapply(X=data$X_train,FUN=dim)
#' sapply(X=data$y_test,FUN=length)
#' sapply(X=data$X_test,FUN=dim)
#' dim(data$beta)
#' 
sim.data.transfer <- function(prob.common=0.05,prob.separate=0.05,q=3,n0=c(50,100,200),n1=10000,p=200,rho=0.5,family="gaussian"){
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

#' @title Data simulation for multi-task learning
#' 
#' @description
#' Simulates data for multi-task learning.
#' 
#' @export
#' 
#' @inheritParams sim.data.transfer
#' 
#' @returns
#' Returns list with slots
#' y_train (n0 x q matrix),
#' X_train (n0 x p matrix),
#' y_test (n1 x q matrix),
#' X_test (n1 x p matrix),
#' and beta (p x q matrix).
#' 
#' @examples
#' data <- sim.data.multiple()
#' sapply(X=data,FUN=dim)
#'  
sim.data.multiple <- function(prob.common=0.05,prob.separate=0.05,q=3,n0=100,n1=10000,p=200,rho=0.5,family="gaussian"){
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

#' @title Calculate deviance
#' 
#' @description
#' Calculates Gaussian deviance (mean-squared error) and binomial deviance.
#' 
#' @export
#' 
#' @param y response
#' @param y_hat predictor
#' @param family character "gaussian" or "binomial"
#' 
#' @examples
#' 
#' n <- 100
#' family <- "gaussian"
#' y <- stats::rnorm(n=n)
#' y_hat <- stats::rnorm(n=n)
#' calc.metric(y=y,y_hat=y_hat,family=family)
#' 
#' family <- "binomial"
#' y <- stats::rbinom(n=n,size=1,prob=0.5)
#' y_hat <- stats::runif(n=n)
#' calc.metric(y=y,y_hat=y_hat,family=family)
#' 
calc.metric <- function(y,y_hat,family){
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

#' @title  Folds for multi-task learning
#' 
#' @export
#' 
#' @param y matrix with n rows (samples) and q columns (outcomes)
#' @param family character "gaussian" or "binomial"
#' @param nfolds integer between 2 and n
#' 
#' @examples
#' 
#' family <- "binomial"
#' y <- sim.data.multiple(family=family)$y_train
#' fold <- make.folds.multi(y=y,family=family)
#' table(fold)
#' table(y[,3],fold)
#' 
make.folds.multi <- function(y,family,nfolds=10){
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
      #warning("Improve this function!")
      if(sum(cond)==1){ # Should this be sum(cond)<=nfolds?
        foldid[cond] <- cands
      } else {
        foldid[cond] <- sample(x=rep(x=cands,length.out=sum(cond))) # Should this be sample(cands)?
      }
    }
  }
  return(foldid)
}

#' @title Folds for transfer learning
#' 
#' @export
#' 
#' @param y list of q numeric vectors
#' @inheritParams make.folds.multi
#' 
#' @examples
#' family <- "binomial"
#' y <- sim.data.transfer(family=family)$y_train
#' fold <- make.folds.trans(y,family=family)
#' 
make.folds.trans <- function(y,family,nfolds=10){
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

get.info <- function(x,y,family){
  if(length(x)!=length(y)){stop("different q")}
  if(any(sapply(X=x,FUN=base::nrow)!=sapply(X=y,FUN=base::length))){stop("different n")}
  if(any(diff(sapply(X=x,FUN=base::ncol))!=0)){stop("different p")}
  q <- length(x)
  n <- sapply(X=x,FUN=base::nrow)
  p <- ncol(x[[1]])
  list <- list(q=q,n=n,p=p)
  return(list)
}

#' @title Data fusion
#' 
#' @export
#'
#' @param x list of q matrices, with n_1,...,n_q rows and with p columns
#' @param y list of q vectors, of length n_1,...,n_q, or NULL (default)
#' @param foldid list of q vectors, of length n_1,...n_q, or NULL (default)
#' 
#' @examples
#' data <- sim.data.transfer()
#' sapply(X=data$y_train,FUN=length)
#' sapply(X=data$X_train,FUN=dim)
#' fuse <- fuse.data(x=data$X_train,y=data$y_train)
#' length(fuse$y)
#' dim(fuse$x)
#' table(fuse$index)
#' 
fuse.data <- function(x,y=NULL,foldid=NULL){
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

# glm.shrink <- function(x,y,family,trial=TRUE){
#   info <- get.info(x=x,y=y,family=family)
#   
#   nfolds <- 10
#   foldid <- make.folds.trans(y=y,family=family,nfolds=nfolds)
#   fuse <- fuse.data(x=x,y=y,foldid=foldid)
#   
#   if(trial){
#     penalty.factor <- 1/abs(stats::cor(x=fuse$x,y=fuse$y,method="spearman"))
#     penalty.factor[is.na(penalty.factor)] <- Inf
#   } else {
#     penalty.factor <- rep(x=1,times=ncol(fuse$x))
#   }
#   glm.com.ext <- glmnet::cv.glmnet(x=fuse$x,y=fuse$y,family=family,penalty.factor=penalty.factor)
#   # To be used above: The adaptive lasso has been implemented in another file.
#   coef.com.ext <- coef(object=glm.com.ext,s="lambda.min")[-1]
#   glm.sep.ext <- list()
#   for(i in seq_len(info$q)){
#     offset <- x[[i]] %*% coef.com.ext
#     glm.sep.ext[[i]] <- glmnet::glmnet(x=x[[i]],y=y[[i]],family=family,offset=offset)
#   }
#   
#   y_hat <- list()
#   for(i in seq_len(info$q)){
#     nlambda <- length(glm.sep.ext[[i]]$lambda)
#     y_hat[[i]] <- matrix(data=NA,nrow=info$n[i],ncol=nlambda)
#   }
#   
#   for(k in seq_len(nfolds)){
#     cond <- fuse$foldid==k
#     if(trial){
#       penalty.factor <- 1/abs(stats::cor(x=fuse$x[!cond,],y=fuse$y[!cond],method="spearman"))
#       penalty.factor[is.na(penalty.factor)] <- Inf
#     } else {
#       penalty.factor <- rep(x=1,times=ncol(fuse$x))
#     }
#     glm.com.int <- glmnet::glmnet(x=fuse$x[!cond,],y=fuse$y[!cond],family=family,penalty.factor=penalty.factor)
#     coef.com.int <- stats::coef(object=glm.com.int,s=glm.com.ext$lambda.min)[-1]
#     glm.sep.int <- list()
#     for(i in seq_len(info$q)){
#       cond <- foldid[[i]]==k
#       offset <- x[[i]] %*% coef.com.int
#       glm.sep.int[[i]] <- glmnet::glmnet(x=x[[i]][!cond,],y=y[[i]][!cond],family=family,offset=offset[!cond])
#       y_hat[[i]][cond,] <- stats::predict(object=glm.sep.int[[i]],newx=x[[i]][cond,],newoffset=offset[cond],s=glm.sep.ext[[i]]$lambda,type="response")
#     }
#   }
#   
#   metric <- list()
#   id.min <- lambda.min <- rep(x=NA,times=info$q)
#   for(i in seq_len(info$q)){
#     metric[[i]] <- apply(X=y_hat[[i]],MARGIN=2,FUN=function(x) calc.metric(y=y[[i]],y_hat=x,family=family))
#     id.min[i] <- which.min(metric[[i]])
#     lambda.min[i] <- glm.sep.ext[[i]]$lambda[id.min[i]]
#   }
#   id.min <- sapply(X=metric,FUN=base::which.min)
#   
#   list <- list(glm.com=glm.com.ext,glm.sep=glm.sep.ext,lambda.min=lambda.min,info=info)
#   class(list) <- "glm.shrink"
#   return(list)  
# }
# 
# predict.glm.shrink <- function(object,newx){
#   q <- length(newx)
#   coef.com <- stats::coef(object=object$glm.com,s="lambda.min")[-1]
#   y_hat <- list()
#   for(i in seq_len(q)){
#     newoffset <- newx[[i]] %*% coef.com
#     y_hat[[i]] <- stats::predict(object=object$glm.sep[[i]],newx=newx[[i]],newoffset=newoffset,s=object$lambda.min[i],type="response")
#   }
#   return(y_hat)
# }
# 
# coef.glm.shrink <- function(object){
#   coef.com <- as.numeric(stats::coef(object=object$glm.com,s="lambda.min"))
#   coef.tot <- matrix(data=NA,nrow=object$info$p+1,ncol=object$info$q)
#   for(i in seq_len(object$info$q)){
#     coef.tot[,i] <- as.numeric(stats::coef(object$glm.sep[[i]],s=object$lambda.min[i])) + coef.com
#   }
#   alpha <- coef.tot[1,]
#   beta <- coef.tot[-1,]
#   list <- list(alpha=alpha,beta=beta)
#   return(list)
# }
# 
# #object <- glm.shrink(x=X_train,y=y_train,family=family)
# 
# 
# # This version uses the same set of weights for each problem.
# combcoef <- function(coef,weight){
#   weight <- weight/sum(weight)
#   #weight <- 0*weight + 1
#   #positive <- apply(X=coef,MARGIN=1,FUN=function(x) sum(weight*pmax(0,x))) # original
#   #negative <- abs(apply(X=coef,MARGIN=1,FUN=function(x) sum(weight*pmin(0,x)))) # original
#   positive <- apply(X=coef,MARGIN=1,FUN=function(x) max(c(0,weight*x)))
#   negative <- apply(X=coef,MARGIN=1,FUN=function(x) abs(min(c(0,weight*x))))
#   #warning("start test")
#   #positive <- apply(X=coef,MARGIN=1,FUN=function(x) max(c(0,x)))
#   #negative <- apply(X=coef,MARGIN=1,FUN=function(x) abs(min(c(0,x))))
#   #warning("end test")
#   lower.limits <- rep(c(0,-Inf),each=nrow(coef))
#   upper.limits <- rep(c(Inf,0),each=nrow(coef))
#   weight <- c(positive,negative)
#   penalty.factor <- 1/weight
#   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight=weight,penalty.factor=penalty.factor)
# }
# 
# 
# comb_split <- function(coef,weight,id){
#   weight <- weight/sum(weight)
#   #--- source ---
#   positive <- apply(X=coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) max(c(0,weight[-id]*x)))
#   negative <- apply(X=coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) abs(min(c(0,weight[-id]*x))))
#   lower.limits <- rep(c(0,-Inf),each=nrow(coef))
#   upper.limits <- rep(c(Inf,0),each=nrow(coef))
#   weight.source <- c(positive,negative)
#   #--- target ---
#   positive <- pmax(0,coef[,id])
#   negative <- abs(pmin(0,coef[,id]))
#   weight.target <- c(positive,negative)
#   
#   #warning("start temporary")
#   if(any(weight.source!=0)){
#     weight.source <- weight.source/sum(weight.source)
#   }
#   if(any(weight.target!=0)){
#     weight.target <- weight.target/sum(weight.target)
#   }
#   #warning("end temporary")
#   
#   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight.source=weight.source,weight.target=weight.target)
# }


comb_split_trial <- function(coef,id){

  # #--- prediction ---
  # metric <- numeric()
  # for(i in seq_len(ncol(coef))){
  #    eta <- x_sta[[id]] %*% coef[,i]
  #    y_hat <- mean_function(eta=eta,family=family)
  #    metric[i] <- calc.metric(y=y_sta[[id]],y_hat=y_hat,family=family)
  # }
  # mean <- rep(x=y_sta[[i]],times=length(y_sta[[i]]))
  # metric.mean <- calc.metric(y=y_sta[[id]],y_hat=y_hat,family=family)
  # cor <- 1*(metric[-id]<metric.mean)

  # #--- correlation ---
  #cor <- stats::cor(x=coef[,id],y=coef[,-id])
  #cor[is.na(cor)] <- 0
  #cor[cor>0] <- 1
  #cor[cor<0] <- 0
  #cor <- matrix(data=cor,nrow=nrow(coef),ncol=length(cor),byrow=TRUE)

  #--- equality ---
  cor <- 1

  #warning("remove binarisation!")
  #coef <- sign(coef) # temporary trial

  #--- external weights ---

  warning("put back max/min (instead of sum)")
  #positive <- apply(X=cor*coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) max(c(0,x)))
  #negative <- apply(X=cor*coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) abs(min(c(0,x))))

  positive <- apply(X=cor*coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) sum(pmax(0,x)))
  negative <- apply(X=cor*coef[,-id,drop=FALSE],MARGIN=1,FUN=function(x) sum(abs(pmin(0,x))))

  lower.limits <- rep(c(0,-Inf),each=nrow(coef))
  upper.limits <- rep(c(Inf,0),each=nrow(coef))
  weight.ext <- c(positive,negative)

  #--- internal weights ---
  positive <- pmax(0,coef[,id])
  negative <- abs(pmin(0,coef[,id]))
  weight.int <- c(positive,negative)

  warning("removed weight standardisation")
  #warning("reactivated weight standardisation")
  #if(any(weight.ext!=0)){
  #  weight.ext <- weight.ext/sum(weight.ext)
  #}
  #if(any(weight.int!=0)){
  #  weight.int <- weight.int/sum(weight.int)
  #}
  #warning("end temporary")

  list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight.source=weight.ext,weight.target=weight.int)
}


# combcoef_shrink <- function(coef,weight){
#   q <- ncol(coef)
#   source.weight <- weight/sum(weight)
#   #cor <- stats::cor(x=coef,method="spearman")
#   pvalue <- matrix(data=0,nrow=q,ncol=q)
#   for(i in 1:(ncol(coef)-1)){
#     for(j in (i+1):ncol(coef)){
#       pvalue[i,j] <- pvalue[j,i] <- stats::cor.test(x=coef[,i],y=coef[,j],method="spearman")$p.value
#     }
#   }
#   signif <- 1*(pvalue<=0.05) # rename to sign (here and below)
#   weight <- penalty.factor <- list()
#   for(i in seq_len(q)){
#     #positive <- apply(X=coef,MARGIN=1,FUN=function(x) sum(source.weight*pmax(0,signif[i,])*pmax(0,x)))
#     #negative <- abs(apply(X=coef,MARGIN=1,FUN=function(x) sum(source.weight*pmax(0,signif[i,])*pmin(0,x))))
#     positive <- apply(X=coef,MARGIN=1,FUN=function(x) max(c(0,signif[i,]*source.weight*x)))
#     negative <- apply(X=coef,MARGIN=1,FUN=function(x) abs(min(c(0,signif[i,]*source.weight*x))))
#     weight[[i]] <- c(positive,negative)
#     penalty.factor[[i]] <- 1/weight[[i]]
#   }
#   lower.limits <- rep(c(0,-Inf),each=nrow(coef))
#   upper.limits <- rep(c(Inf,0),each=nrow(coef))
#   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight=weight,penalty.factor=penalty.factor)
# }
# # 
# # combcoef_corshrink <- function(coef,weight){
# #   q <- ncol(coef)
# #   shrink <- t(apply(coef,1,function(x) CorShrink::CorShrinkVector(corvec=x,nsamp_vec=weight)))
# #   
# #   weight <- weight/sum(weight)
# #   
# #   penalty.factor <- list()
# #   for(i in seq_len(q)){
# #     positive <- apply(X=coef,MARGIN=1,FUN=function(x) sum(weight*pmax(0,cor[i,])*pmax(0,x)))
# #     negative <- abs(apply(X=coef,MARGIN=1,FUN=function(x) sum(weight*pmax(0,cor[i,])*pmin(0,x))))
# #     penalty.factor[[i]] <- 1/c(positive,negative)
# #   }
# #   lower.limits <- rep(c(0,-Inf),each=nrow(coef))
# #   upper.limits <- rep(c(Inf,0),each=nrow(coef))
# #   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,penalty.factor=penalty.factor)
# # }
# 
# # This version imposes uniform signs.
# combcoef_trial <- function(coef,weight){
#   weight <- weight/sum(weight)
#   positive <- apply(X=coef,MARGIN=1,FUN=function(x) sum(weight*pmax(0,x)))
#   negative <- abs(apply(X=coef,MARGIN=1,FUN=function(x) sum(weight*pmin(0,x))))
#   lower.limits <- rep(c(0,-Inf),each=nrow(coef))
#   upper.limits <- rep(c(Inf,0),each=nrow(coef))
#   penalty.factor <- 1/c(positive,negative)
#   penalty.factor[c(negative>positive,positive>negative)] <- 0
#   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,penalty.factor=penalty.factor)
# }
# 
# # This version uses a different set of weights for each problem.
# combcoef2 <- function(alpha,beta,y,x,family){
#   q <- length(x)
#   p <- ncol(x[[1]])
#   n <- sapply(X=y,FUN=base::length)
#   n <- n/sum(n)
#   metric <- matrix(data=NA,nrow=q,ncol=q+1)
#   for(i in seq_len(q)){
#     for(j in seq_len(q)){
#       eta_hat <- alpha[j] + x[[i]] %*% beta[,j]
#       y_hat <- mean_function(eta=eta_hat,family=family)
#       metric[i,j] <- calc.metric(y=y[[i]],y_hat=y_hat,family=family)
#     }
#     metric[i,j+1] <- calc.metric(y=y[[i]],y_hat=rep(x=mean(y[[i]]),times=n[i]),family=family)
#   }
#   metric <- metric[,-(q+1)]/metric[,q+1]
#   select <- apply(X=metric,MARGIN=2,FUN=function(x) 1*(x<1))
#   penalty.factor <- weights <- list()
#   for(i in seq_len(q)){
#     weight <- select[i,]
#     positive <- apply(X=beta,MARGIN=1,FUN=function(x) max(0,weight*n*x))
#     negative <- abs(apply(X=beta,MARGIN=1,FUN=function(x) min(0,weight*n*x)))
#     weights[[i]] <- c(positive,negative)
#     penalty.factor[[i]] <- 1/weight[[i]]
#   }
#   lower.limits <- rep(c(0,-Inf),each=nrow(beta))
#   upper.limits <- rep(c(Inf,0),each=nrow(beta))
#   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight=weights,penalty.factor=penalty.factor)
#   return(list)
# }
# 
# combpval <- function(pval.inc,pval.dec){
#   p <- nrow(pval.inc)
#   q <- ncol(pval.dec)
#   #prob <- apply(rbind(pval.inc,pval.dec),2,function(x) mean(x<=0.05))
#   #pval.inc[] <- 
#   #pval.dec[] <- 
#   comb.inc <- apply(pval.inc,1,function(x) palasso:::.combine(x,method="fisher"))
#   comb.dec <- apply(pval.dec,1,function(x) palasso:::.combine(x,method="fisher"))
#   pvalue <- c(comb.inc,comb.dec)
#   weight <- -log10(pvalue)
#   weight[weight<=1e-06] <- 0
#   # CONSIDER CROSS-VALIDATING EXPONENTS!
#   penalty.factor <- 1/weight
#   penalty.factor[weight<=1e-06] <- Inf
#   lower.limits <- rep(c(0,-Inf),each=p)
#   upper.limits <- rep(c(Inf,0),each=p)
#   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight=weight,penalty.factor=penalty.factor)
#   return(list)
# }
# 
# signpval <- function(x,y){
#   # This function is under development.
#   temp <- apply(x,2,function(x) stats::cor.test(x,y,method="spearman",exact=FALSE))
#   temp.sign <- sapply(temp,function(x) sign(x$estimate))
#   temp.pval <- -log10(sapply(temp,function(x) x$p.value))
#   prior <- temp.sign*temp.pval
#   return(prior)
# }
# 
# onesidecor <- function(x,y){
#   inc <- pmin(1,2*apply(x,MARGIN=2,FUN=function(x) stats::cor.test(x=x,y=y,method="spearman",exact=FALSE,alternative="greater")$p.value))
#   dec <- pmin(1,2*apply(X=x,MARGIN=2,FUN=function(x) stats::cor.test(x=x,y=y,method="spearman",exact=FALSE,alternative="less")$p.value))
#   #warning("start temporary")
#   #inc[inc>0.05] <- 1 # trial 2024-07-10
#   #dec[dec>0.05] <- 1 # trial 2024-07-10
#   #warning("end temporary")
#   list <- list(inc=inc,dec=dec)
#   return(list)
# }
# 
# comb_calib <- function(y,x,coef,n,family){
#   n <- n/sum(n)
#   q <- length(x)
#   p <- ncol(x[[1]])
#   # calibration
#   prior <- list()
#   for(i in seq_len(q)){
#     prior[[i]] <- matrix(data=NA,nrow=p,ncol=q)
#     for(j in seq_len(q)){
#       prior[[i]][,j] <- transreg:::.exp.multiple(y=y[[i]],X=x[[i]],family=family,prior=coef[,j,drop=FALSE])$beta
#     }
#   }
#   # combination
#   weight <- penalty.factor <- list()
#   for(i in seq_len(q)){
#     positive <- apply(X=prior[[i]],MARGIN=1,FUN=function(x) max(c(0,n*x)))
#     negative <- apply(X=prior[[i]],MARGIN=1,FUN=function(x) abs(min(c(0,n*x))))
#     weight[[i]] <- c(positive,negative)
#     penalty.factor[[i]] <- 1/weight[[i]]
#   }
#   lower.limits <- rep(c(0,-Inf),each=p)
#   upper.limits <- rep(c(Inf,0),each=p)
#   list <- list(lower.limits=lower.limits,upper.limits=upper.limits,weight=weight,penalty.factor=penalty.factor)
#   return(list)
# }
# 
# sparselink <- function(x,y,family,alpha.one=0.95,alpha.two=1,drop=FALSE,tune=TRUE){ # Change back to alpha.one=NA (correlation). Under tune=TRUE, alpha.one=Inf (transformed p-values) and alpha.one=0 (ridge) work equally well.
#   
#   if(!is.na(alpha.one)){warning("Not using correlation!")}
#   
#   #drop <- FALSE # was test <- FALSE # Use this to check whether a source is used.
#   
#   n <- sapply(X=y,FUN=base::length)
#   p <- ncol(x[[1]])
#   q <- length(x)
#   
#   nfolds <- 10
#   foldid <- make.folds.trans(y=y,family=family,nfolds=nfolds)
#   
#   # overall fit
#   glm.one.ext <- glm.two.ext <- list()
#   alpha.ext <- rep(x=NA,times=q)
#   beta.ext <- pval.ext.inc <- pval.ext.dec <- matrix(data=NA,nrow=p,ncol=q)
#   for(i in seq_len(q)){
#     if(is.na(alpha.one)){
#       beta.ext[,i] <- stats::cor(x=x[[i]],y=y[[i]],method="spearman")
#       beta.ext[,i][is.na(beta.ext[,i])] <- 0
#       #warning("start temporary")
#       #beta.ext[,i] <- signpval(x=x[[i]],y=y[[i]])
#       #warning("end temporary")
#     } else if(is.infinite(alpha.one)){
#       temp <- onesidecor(x=x[[i]],y=y[[i]])
#       pval.ext.inc[,i] <- temp$inc
#       pval.ext.dec[,i] <- temp$dec
#     } else {
#       glm.one.ext[[i]] <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],family=family,alpha=alpha.one) # use fixed foldid?
#       temp <- coef(object=glm.one.ext[[i]],s="lambda.min")
#       alpha.ext[i] <- temp[1]
#       beta.ext[,i] <- temp[-1]
#     }
#   }
#   
#   if(is.infinite(alpha.one)){
#     rest.ext <- combpval(pval.inc=pval.ext.inc,pval.dec=pval.ext.dec)
#   } else if(drop){
#     rest.ext <- combcoef2(alpha=alpha.ext,beta=beta.ext,y=y,x=x,family=family)
#     #rest.ext <- combcoef_shrink(coef=beta.ext,weight=n)
#     #rest.ext <- comb_calib(y=y,x=x,coef=beta.ext,n=n,family=family)
#   } else {
#     rest.ext <- combcoef(coef=beta.ext,weight=n)
#   }
#   
#   exponent <- seq(from=0,to=2,by=0.25) # verify range
#   #exponent <- seq(from=0,to=1,by=0.2) # trial
#   
#   for(i in seq_len(q)){
#     if(tune){
#       if(drop){
#         weight <- rest.ext$weight[[i]]
#       } else {
#         weight <- rest.ext$weight
#       }
#     } else {
#       if(drop){
#         penalty.factor <- rest.ext$penalty.factor[[i]]
#       } else {
#         penalty.factor <- rest.ext$penalty.factor
#       }
#     }
#     if(tune){
#       glm.two.ext[[i]] <- list()
#       for(l in seq_along(exponent)){
#         if(all(weight^exponent[l]==0)){
#           glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lambda=99e99,alpha=alpha.two)
#         } else {
#           glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lower.limits=rest.ext$lower.limits,upper.limits=rest.ext$upper.limits,penalty.factor=1/(weight^exponent[l]),alpha=alpha.two)
#         }
#       }
#     } else {
#       glm.two.ext[[i]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lower.limits=rest.ext$lower.limits,upper.limits=rest.ext$upper.limits,penalty.factor=penalty.factor,alpha=alpha.two)
#     }
#   }
#   
#   # cross-validation
#   y_hat <- list()
#   if(tune){
#     for(i in seq_len(q)){
#       y_hat[[i]] <- list()
#       for(l in seq_along(exponent)){
#         nlambda <- length(glm.two.ext[[i]][[l]]$lambda)
#         y_hat[[i]][[l]] <- matrix(data=NA,nrow=n[i],ncol=nlambda)
#       }
#     }
#   } else {
#     for(i in seq_len(q)){
#       nlambda <- length(glm.two.ext[[i]]$lambda)
#       y_hat[[i]] <- matrix(data=NA,nrow=n[i],ncol=nlambda)
#     }
#   }
#   
#   for(k in seq_len(nfolds)){
#     alpha.int <- rep(x=NA,times=q)
#     beta.int <- matrix(data=NA,nrow=p,ncol=q)
#     pval.int.inc <- pval.int.dec <- matrix(data=NA,nrow=p,ncol=q)
#     for(i in seq_len(q)){
#       cond <- foldid[[i]]==k
#       if(is.na(alpha.one)){
#         beta.int[,i] <- stats::cor(x=x[[i]][!cond,],y=y[[i]][!cond],method="spearman")
#         beta.int[,i][is.na(beta.int[,i])] <- 0
#         #warning("start temporary")
#         #beta.int[,i] <- signpval(x=x[[i]][!cond,],y=y[[i]][!cond])
#         #warning("end temporary")
#       } else if(is.infinite(alpha.one)){
#         temp <- onesidecor(x=x[[i]][!cond,],y=y[[i]][!cond])
#         pval.int.inc[,i] <- temp$inc
#         pval.int.dec[,i] <- temp$dec
#       } else {
#         glm.one.int <- glmnet::glmnet(x=x[[i]][!cond,],y=y[[i]][!cond],family=family,alpha=alpha.one)
#         temp <- coef(object=glm.one.int,s=glm.one.ext[[i]]$lambda.min)
#         alpha.int[i] <- temp[1]
#         beta.int[,i] <- temp[-1]
#       }
#     }
#     
#     if(is.infinite(alpha.one)){
#       rest.int <- combpval(pval.inc=pval.int.inc,pval.dec=pval.int.dec)
#     } else if(drop){
#       y_temp <- x_temp <- list()
#       for(i in seq_len(q)){
#         y_temp[[i]] <- y[[i]][foldid[[i]]==k]
#         x_temp[[i]] <- x[[i]][foldid[[i]]==k,]
#       }
#       rest.int <- combcoef2(alpha=alpha.int,beta=beta.int,y=y_temp,x=x_temp,family=family)
#       #rest.int <- combcoef_shrink(coef=beta.int,weight=n) 
#       #rest.int <- comb_calib(y=y,x=x,coef=beta.int,n=n,family=family)
#     } else {
#       rest.int <- combcoef(coef=beta.int,weight=n)
#     }
#     
#     for(i in seq_len(q)){
#       cond <- foldid[[i]]==k
#       if(tune){
#         if(drop){
#           weight <- rest.int$weight[[i]]
#         } else {
#           weight <- rest.int$weight
#         }
#       } else {
#         if(drop){
#           penalty.factor <- rest.int$penalty.factor[[i]]
#         } else {
#           penalty.factor <- rest.int$penalty.factor
#         }
#       }
#       if(tune){
#         for(l in seq_along(exponent)){
#           if(all(weight^exponent[l]==0)){
#             # This should be handled differently.
#             #y_hat[[i]][[l]][cond,] <- predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=Inf,type="response")
#             glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lambda=99e99,alpha=alpha.two)
#           } else {
#             glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lower.limits=rest.int$lower.limits,upper.limits=rest.int$upper.limits,penalty.factor=1/(weight^exponent[l]),alpha=alpha.two)
#           }
#           y_hat[[i]][[l]][cond,] <- predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=glm.two.ext[[i]][[l]]$lambda,type="response")
#         }
#       } else {
#         glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lower.limits=rest.int$lower.limits,upper.limits=rest.int$upper.limits,penalty.factor=penalty.factor,alpha=alpha.two)
#         y_hat[[i]][cond,] <- predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=glm.two.ext[[i]]$lambda,type="response")
#       }
#       
#     }
#   }
#   
#   if(tune){
#     metric <- list()
#     id.min <- lambda.min <- rep(x=NA,times=q)
#     for(i in seq_len(q)){
#       metric[[i]] <- list()
#       for(l in seq_along(exponent)){
#         metric[[i]][[l]] <- apply(X=y_hat[[i]][[l]],MARGIN=2,FUN=function(x)
#           calc.metric(y=y[[i]],y_hat=x,family=family))
#       }
#       min <- sapply(metric[[i]],min)
#       tryCatch(expr=graphics::plot(x=exponent,y=min,type="o"),error=function(x) NULL)
#       id.exp <- which.min(min)
#       id.min[i] <- which.min(metric[[i]][[id.exp]])
#       lambda.min[i] <- glm.two.ext[[i]][[id.exp]]$lambda[id.min[i]]
#       glm.two.ext[[i]] <- glm.two.ext[[i]][[id.exp]]
#     }
#     
#   } else {
#     metric <- list()
#     id.min <- lambda.min <- rep(x=NA,times=q)
#     for(i in seq_len(q)){
#       metric[[i]] <- apply(X=y_hat[[i]],MARGIN=2,FUN=function(x)
#         calc.metric(y=y[[i]],y_hat=x,family=family))
#       id.min[i] <- which.min(metric[[i]])
#       lambda.min[i] <- glm.two.ext[[i]]$lambda[id.min[i]]
#     }
#   }
#   
#   list <- list(glm.one=glm.one.ext,glm.two=glm.two.ext,lambda.min=lambda.min,info=data.frame(p=p,q=q))
#   class(list) <- "sparselink"
#   return(list)
# }

#' @export
predict.sparselink <- function(object,newx,weight=NULL){
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

#' @export
coef.sparselink <- function(object){
  id <- object$weight.ind
  p <- object$info$p
  q <- object$info$q #length(object$lambda.min)
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


# #object <- sparselink(x=X_train,y=y_train,family=family)
# #y_hat <- predict(object=object,newx=X_test)
# #coef <- coef(object=object)

glm.common <- function(x,y,family,alpha=1){
  family <- unique(family)
  if(length(family)>1){stop("requires unique family")}
  object <- list()
  object$info <- get.info(x=x,y=y,family=family)
  fuse <- fuse.data(x=x,y=y,foldid=NULL)
  object$cv.glmnet <- glmnet::cv.glmnet(x=fuse$x,y=fuse$y,family=family,alpha=alpha)
  class(object) <- "glm.common"
  return(object)
}

#' @export
predict.glm.common <- function(object,newx){
  fuse <- fuse.data(x=newx,y=NULL,foldid=NULL)
  temp <- stats::predict(object=object$cv.glmnet,newx=fuse$x,s=object$cv.glmnet$lambda.min,type="response")
  y_hat <- tapply(X=temp,INDEX=fuse$index,FUN=function(x) x)
  return(y_hat)
}

#' @export
coef.glm.common <- function(object){
  coef <- stats::coef(object=object$cv.glmnet,s="lambda.min")
  alpha <- rep(x=coef[1],times=object$info$q)
  beta <- matrix(data=coef[-1],nrow=object$info$p,ncol=object$info$q)
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

#object <- glm.common(x=X_train,y=y_train,family=family)
#y_hat <- predict(object,newx=X_test)
#coef <- coef(object)

glm.separate <- function(x,y,family,alpha=1){
  if(is.matrix(x) & is.matrix(y)){
    q <- ncol(y)
    x <- replicate(n=q,expr=x,simplify=FALSE)
    y <- apply(y,2,function(x) x,simplify=FALSE)
  }
  if(length(family)==1){
    family <- rep(x=family,times=length(y))
  }
  object <- list()
  object$info <- get.info(x=x,y=y,family=family)
  object$cv.glmnet <- list()
  for(i in seq_len(object$info$q)){
    object$cv.glmnet[[i]] <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],family=family[i],alpha=alpha)
  }
  class(object) <- "glm.separate"
  return(object)
}

#' @export
predict.glm.separate <- function(object,newx){
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

#' @export
coef.glm.separate <- function(object){
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

#object <- glm.separate(x=X_train,y=y_train,family=family)
#y_hat <- predict(object,newx=X_test)
#coef <- coef(object)

if(FALSE){
  # family="mgaussian"
  n <- 100
  p <- 200
  q <- 3
  x <- matrix(data=stats::rnorm(n=n*p),nrow=n,ncol=p)
  theta <- stats::rbinom(n=p,size=1,prob=0.05)*stats::rnorm(n=p)
  beta <- matrix(data=NA,nrow=p,ncol=q)
  for(i in seq_len(q)){
    beta[,i] <- theta + stats::rbinom(n=p,size=1,prob=0.05)*stats::rnorm(n=p)
  }
  y <- x %*% beta + stats::rnorm(n*q)
  object <- glmnet::cv.glmnet(x=x,y=y,family="mgaussian")
  sapply(coef(object,s="lambda.min"),function(x) as.numeric(x))
  
  object <- glm.mgaussian(x=x,y=y,family="gaussian")
  predict(object,newx=x)
}

glm.mgaussian <- function(x,y,family,alpha){
  object <- list()
  family <- unique(family)
  if(length(family)>1){stop("requires unique family")}
  if(family!="gaussian"){stop("requires gaussian")}
  object$cv.glmnet <- glmnet::cv.glmnet(x=x,y=y,family="mgaussian",alpha=alpha)
  class(object) <- "glm.mgaussian"
  return(object)
}

#' @export
predict.glm.mgaussian <- function(object,newx){
  y_hat <- stats::predict(object$cv.glmnet,newx=newx,s="lambda.min")
  apply(y_hat,2,function(x) x,simplify=FALSE)
}

#' @export
coef.glm.mgaussian <- function(object){
  coef <- stats::coef(object$cv.glmnet,s="lambda.min")
  alpha <- sapply(coef,function(x) x[1])
  beta <- sapply(coef,function(x) x[-1])
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

glm.transfer <- function(x,y,family,alpha=1){
  family <- unique(family)
  if(length(family)>1){stop("glmtrans requires unique family")}
  q <- length(x)
  source <- list()
  for(i in seq_len(q)){
    source[[i]] <- list(y=y[[i]],x=x[[i]])
  }
  object <- list()
  for(i in seq_len(q)){
    invisible(utils::capture.output(object[[i]] <- glmtrans::glmtrans(target=list(y=y[[i]],x=x[[i]]),source=source[-i],family=family,alpha=alpha)))
  }
  class(object) <- "glm.transfer"
  return(object)
}

#' @export
predict.glm.transfer <- function(object,newx){
  q <- length(newx)
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- stats::predict(object=object[[i]],newx=newx[[i]],type="response",s="lambda.min")
  }
  return(y_hat)
}

#' @export
coef.glm.transfer <- function(object){
  coef <- sapply(object,function(x) x$beta)
  alpha <- coef[1,]
  beta <- coef[-1,]
  list <- list(alpha=alpha,beta=beta)
  return(list)
}

#object <- glm.transfer(x=X_train,y=y_train,family=family)
#y_hat <- predict(object,newx=X_train)
#coef(object)

#cor(y_train[[3]],y_hat[[3]])

#' @title Train and test model
#'
#' @description
#' Trains and test prediction models
#' 
#' @param y_train target of training samples: vector of length n
#' @param X_train features of training samples: n x p matrix
#' @param y_test target of testing samples: vector of length m
#' @param X_test features of testing samples m x p matrix
#' @param family character "gaussian" or "binomial"
#' @param alpha.init elastic net mixing parameter for initial regressions
#' @param alpha elastic net mixing paramter
#' @param method character vector
#' @param type character
#'
traintest <- function(y_train,X_train,y_test=NULL,X_test=NULL,family,alpha,method=c("glm.separate","glm.transfer","sparselink","glm.common"),alpha.init,type){
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
    #func <- eval(parse(text=paste0("glm.",method[i])))
    func <- eval(parse(text=paste0(method[i])))
    start <- Sys.time()
    hyperpar <- NULL
    if(method[i]=="sparselink"){
      object <- func(x=X_train,y=y_train,family=family,alpha.init=alpha.init,alpha=alpha,type=type)
      hyperpar <- object$weight.min
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
        deviance[i,j] <- calc.metric(y=y_test[[j]],y_hat=y_hat[[i]][[j]],family=family[j])
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


# This cross-validation function only works for transfer learning (not for multi-task learning).
crossval <- function(y,X,family,alpha=1,nfolds=10,method=c("glm.separate","glm.transfer","sparselink","glm.common"),alpha.init,type){
  #if(is.matrix(y) & is.matrix(X)){
  #  mode <- "multiple"
  #  foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
  #  n <- nrow(y)
  #  q <- ncol(y)
  #  y <- apply(y,2,function(x) x,simplify=FALSE)
  #  X <- replicate(n=q,expr=X,simplify=FALSE)
  #  foldid <- replicate(n=q,expr=foldid,simplify=FALSE)
  #} else {
    mode <- "transfer"
    foldid <- make.folds.trans(y=y,family=family,nfolds=nfolds)
    n <- length(y[[1]])
    q <- length(y)
  #} 
  
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
    test <- traintest(y_train=y_train,X_train=X_train,X_test=X_test,family=family,method=method,alpha=alpha,alpha.init=alpha.init,type=type)
    for(j in seq_len(q)){
      for(k in method){
        y_hat[[j]][foldid[[j]]==i,k] <- test$y_hat[[k]][[j]]
      }
    }
  }
  deviance <- auc <- matrix(data=NA,nrow=length(y),ncol=length(method),dimnames=list(names(y),method))
  for(j in seq_along(y)){
    for(k in method){
      deviance[j,k] <- calc.metric(y=y[[j]],y_hat=y_hat[[j]][,k],family=family)
      if(family=="binomial"){
        auc[j,k] <- pROC::auc(response=y[[j]],predictor=as.vector(y_hat[[j]][,k]),direction="<",levels=c(0,1))
      }
    }
  }
  # refit model on all folds
  cat("refit on all folds","\n")
  refit <- traintest(y_train=y,X_train=X,family=family,method=method,alpha.init=alpha.init,type=type,alpha=alpha)
  list <- list(deviance=deviance,auc=auc,refit=refit)
  return(list)
}



# This cross-validation function only works for transfer learning (not for multi-task learning).
cv.multiple <- function(y,X,family,alpha=1,nfolds=10,method=c("glm.separate","glm.mgaussian","sparselink"),alpha.init,type){
  mode <- "multiple"
  foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
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
    test <- traintest(y_train=y_train,X_train=X_train,X_test=X_test,y_test=y_test,family=family,method=method,alpha=alpha,alpha.init=alpha.init,type=type)
    for(k in method){
        y_hat[[k]][foldid==i,] <- as.matrix(as.data.frame(test$y_hat[[k]]))
    }
  }
  
  deviance <- auc <- matrix(data=NA,nrow=ncol(y),ncol=length(method),dimnames=list(names(y),method))
  for(j in seq_len(ncol(y))){
    for(k in method){
      deviance[j,k] <- calc.metric(y=y[,j],y_hat=y_hat[[k]][,j],family=family)
      if(family=="binomial"){
        auc[j,k] <- pROC::auc(response=y[,j],predictor=as.vector(y_hat[[k]][,j]),direction="<",levels=c(0,1))
      }
    }
  }
  # refit model on all folds
  #cat("refit on all folds","\n")
  #refit <- traintest(y_train=y,X_train=X,family=family,method=method,alpha.init=alpha.init,type=type,alpha=alpha)
  list <- list(deviance=deviance)
  return(list)
}



#test <- traintest(y_train,X_train,y_test,X_test,family)
#test <- crossval(y=y_train,X=X_train,family=family)

#truth <- sample(x=c(-1,0,1),size=20,replace=TRUE)
#estim <- sample(x=c(-1,0,1),size=20,replace=TRUE)
#table(truth,estim)

#' @title Metrics for sign detection
#' 
#' @param truth n times p matrix with entries in -1, 0, 1
#' @param estim n times p matrix with entries in -1, 0, 1
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

#' @title Metrics for sign detection
#' 
#' @param truth vector of length p with entries in -1, 0, 1
#' @param estim vector of length p with entries -1, 0, 1
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

#' @title Plot pairwise differences
#' 
#' @param x setting: character vector
#' @param y0 old value: numeric vector
#' @param y1 new value: numeric vector
#' @param main title
#' @param increase change to arrow NULL, up, down
#' @param cex.axis numeric
#' @param cex.main numeric
#' 
change <- function(x,y0,y1,main="",cex.axis=0.5,cex.main=1,increase=TRUE){
  unique <- unique(x)
  #graphics::par(mfrow=c(1,1),mar=c(3,3,1,1))
  graphics::plot.new()
  xlim <- c(1-0.2,length(unique)+0.2)
  ylim <- range(c(y0,y1),na.rm=TRUE)
  graphics::plot.window(xlim=xlim,ylim=ylim)
  #if(is.list(unique)){
  #  for(i in seq_along(unique)){
  #    #graphics::axis(side=1,at=i,labels=unique[[i]],tick=FALSE,line=0.5,cex.axis=cex.axis)
  #    graphics::mtext(text=unique[[i]][[1]],side=1,line=0.2,at=i,cex=cex.axis)
  #    graphics::mtext(text=unique[[i]][[2]],side=1,line=1.2,at=i,cex=cex.axis)
  #  }
  #} else {
    #graphics::axis(side=1,at=seq_along(unique),labels=unique,tick=FALSE,line=0,cex.axis=cex.axis)
    graphics::mtext(text=unique,side=1,at=seq_along(unique),line=1,cex=cex.axis)  
  #}
  #graphics::mtext(text="common",at=0.5,side=2,padj=1)
  graphics::axis(side=2,cex.axis=cex.axis)
  for(i in seq_along(unique)){
    #if(is.list(unique)){
    #  cond <- sapply(X=x,FUN=function(x) x[[1]]==unique[[i]][[1]] & x[[2]]==unique[[i]][[2]])
    #} else {
      cond <- x==unique[i]
    #}
    graphics::segments(x0=i-0.1,y0=y0[cond],x1=i+0.1,y1=y1[cond],col="grey")
    graphics::points(x=rep(i-0.1,times=sum(cond)),y=y0[cond],col="red",pch=16,cex=0.8)
    graphics::points(x=rep(i+0.1,times=sum(cond)),y=y1[cond],col="blue",pch=16,cex=0.8)
    #pvalue <- stats::wilcox.test(x=y0,y=y1,paired=TRUE)$p.value
  }
  graphics::title(main=main,line=0,cex.main=cex.main)
  #graphics::title(ylab=main,line=3,font=2)
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
  #graphics::arrows(x0=pos,y0=inferior+margin.inferior,
  #                 y1=superior+margin.superior,lwd=2,length=0.08,col="grey")
  #graphics::text(x=pos,y=inferior,labels="-",col="red",font=2,cex=1.2)
  #graphics::text(x=pos,y=superior,labels="+",col="blue",font=2,cex=1.2)
  graphics::arrows(x0=pos,y0=inferior+2*margin.inferior,
                   y1=superior+2*margin.superior,lwd=2,length=0.08,col="grey")
  graphics::text(x=pos,y=inferior+1*margin.inferior,labels="-",col="red",font=2,cex=1.2)
  graphics::text(x=pos,y=superior+1*margin.superior,labels="+",col="blue",font=2,cex=1.2)
  graphics::par(xpd=FALSE)
}

### --- development ---
# 
# glm.tidy <- function(x,y,family,alpha.one=0.95,alpha.two=1){
#   
#   if(FALSE){
#     family <- "binomial"
#     alpha.one <- 0.95
#     alpha.two <- 1
#   }
# 
#   n <- sapply(X=y,FUN=length)
#   p <- ncol(x[[1]])
#   q <- length(x)
# 
#   nfolds <- 10
#   foldid <- make.folds.trans(y=y,family=family,nfolds=nfolds)
# 
#   # overall fit
#   glm.one.ext <- glm.two.ext <- list()
#   alpha.ext <- rep(x=NA,times=q)
#   beta.ext <- matrix(data=NA,nrow=p,ncol=q)
#   for(i in seq_len(q)){
#     glm.one.ext[[i]] <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],family=family,alpha=alpha.one,foldid=foldid[[i]])
#     coef.ext <- coef(object=glm.one.ext[[i]],s="lambda.min")
#     alpha.ext[i] <- coef.ext[1]
#     beta.ext[,i] <- coef.ext[-1]
#   }
#   
#   rest.ext <- combcoef(coef=beta.ext,weight=n)
#   weight <- rest.ext$weight
#   
#   exponent <- seq(from=0,to=2,by=0.25)
#   
#   for(i in seq_len(q)){
#     glm.two.ext[[i]] <- list()
#     for(l in seq_along(exponent)){
#       penalty.factor <- weight^(exponent[l])
#       if(any(is.finite(penalty.factor))){
#         glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lower.limits=rest.ext$lower.limits,upper.limits=rest.ext$upper.limits,penalty.factor=penalty.factor,alpha=alpha.two)
#       }
#       if(all(is.infinite(penalty.factor))||length(glm.two.ext[[i]][[l]]$lambda)==1){
#         glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lambda=99e99,alpha=alpha.two)
#       }
#     }
#   }
# 
#   # cross-validation
#   y_hat <- list()
#   for(i in seq_len(q)){
#     y_hat[[i]] <- list()
#     for(l in seq_along(exponent)){
#       nlambda <- length(glm.two.ext[[i]][[l]]$lambda)
#       y_hat[[i]][[l]] <- matrix(data=NA,nrow=n[i],ncol=nlambda)
#     }
#   }
# 
#   for(k in seq_len(nfolds)){
#     for(i in seq_len(q)){
#       cond <- foldid[[i]]==k
#       
#       alpha.int <- alpha.ext
#       beta.int <- beta.ext
# 
#       glm.one.int <- glmnet::glmnet(x=x[[i]][!cond,],y=y[[i]][!cond],family=family,alpha=alpha.one)
#       coef.int <- coef(object=glm.one.int,s=glm.one.ext[[i]]$lambda.min)
#       alpha.int[i] <- coef.int[1]
#       beta.int[,i] <- coef.int[-1]
#       
#       rest.int <- combcoef(coef=beta.int,weight=n)
#       weight <- rest.int$weight
#     
#       for(l in seq_along(exponent)){
#         penalty.factor <- weight^(exponent[l])
#         if(any(is.finite(penalty.factor))){
#           glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lower.limits=rest.int$lower.limits,upper.limits=rest.int$upper.limits,penalty.factor=penalty.factor,alpha=alpha.two)
#         }
#         if(all(is.infinite(penalty.factor))||length(glm.two.int$lambda)==1){
#           glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lambda=99e99,alpha=alpha.two)
#         }
#         y_hat[[i]][[l]][cond,] <- predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=glm.two.ext[[i]][[l]]$lambda,type="response")
#       }
#     }
#   }
#   
#   metric <- list()
#   id.min <- lambda.min <- rep(x=NA,times=q)
#   for(i in seq_len(q)){
#     metric[[i]] <- list()
#     for(l in seq_along(exponent)){
#       metric[[i]][[l]] <- apply(X=y_hat[[i]][[l]],MARGIN=2,FUN=function(x)
#       calc.metric(y=y[[i]],y_hat=x,family=family))
#     }
#     min <- sapply(metric[[i]],min)
#     tryCatch(expr=graphics::plot(x=exponent,y=min,type="o"),error=function(x) NULL)
#     id.exp <- which.min(min)
#     id.min[i] <- which.min(metric[[i]][[id.exp]])
#     lambda.min[i] <- glm.two.ext[[i]][[id.exp]]$lambda[id.min[i]]
#     glm.two.ext[[i]] <- glm.two.ext[[i]][[id.exp]]
#   }
#     
#   list <- list(glm.one=glm.one.ext,glm.two=glm.two.ext,lambda.min=lambda.min,info=data.frame(p=p,q=q))
#   class(list) <- "sparselink"
#   return(list)
# }

# 
# # This was a working version (2024-07-15).
# glm.retry <- function(x,y,family,alpha.one=0.95,alpha.two=1){
#   
#   n <- sapply(X=y,FUN=base::length)
#   p <- ncol(x[[1]])
#   q <- length(x)
#   
#   nfolds <- 10
#   foldid <- make.folds.trans(y=y,family=family,nfolds=nfolds)
#   
#   # overall fit
#   glm.one.ext <- glm.two.ext <- list()
#   alpha.ext <- rep(x=NA,times=q)
#   beta.ext <- pval.ext.inc <- pval.ext.dec <- matrix(data=NA,nrow=p,ncol=q)
#   for(i in seq_len(q)){
#     glm.one.ext[[i]] <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],family=family,alpha=alpha.one) # use fixed foldid?
#     temp <- coef(object=glm.one.ext[[i]],s="lambda.min")
#     alpha.ext[i] <- temp[1]
#     beta.ext[,i] <- temp[-1]
#   }
#   
#   rest.ext <- combcoef(coef=beta.ext,weight=n)
#   
#   exponent <- seq(from=0,to=2,by=0.25)
#   
#   for(i in seq_len(q)){
#     weight <- rest.ext$weight
#     glm.two.ext[[i]] <- list()
#     for(l in seq_along(exponent)){
#       if(all(weight^exponent[l]==0)){
#         glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lambda=99e99,alpha=alpha.two)
#       } else {
#         glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lower.limits=rest.ext$lower.limits,upper.limits=rest.ext$upper.limits,penalty.factor=1/(weight^exponent[l]),alpha=alpha.two)
#       }
#     }
#   }
#   
#   # cross-validation
#   y_hat <- list()
#   for(i in seq_len(q)){
#     y_hat[[i]] <- list()
#     for(l in seq_along(exponent)){
#       nlambda <- length(glm.two.ext[[i]][[l]]$lambda)
#       y_hat[[i]][[l]] <- matrix(data=NA,nrow=n[i],ncol=nlambda)
#     }
#   }
#   
#   # Changed order of loops on 2024-07-11:
#   
#   #for(k in seq_len(nfolds)){ # original
#   #  alpha.int <- alpha.ext # original
#   #  beta.int <- beta.ext # original
#   #  for(i in seq_len(q)){ # original
#   
#   for(i in seq_len(q)){ # trial   
#     alpha.int <- alpha.ext # trial
#     beta.int <- beta.ext # trial
#     for(k in seq_len(nfolds)){ # trial
#       
#       cond <- foldid[[i]]==k
#       glm.one.int <- glmnet::glmnet(x=x[[i]][!cond,],y=y[[i]][!cond],family=family,alpha=alpha.one)
#       temp <- coef(object=glm.one.int,s=glm.one.ext[[i]]$lambda.min)
#       alpha.int[i] <- temp[1]
#       beta.int[,i] <- temp[-1]
#       
#       rest.int <- combcoef(coef=beta.int,weight=n)
#       
#       cond <- foldid[[i]]==k
#       weight <- rest.int$weight
#       for(l in seq_along(exponent)){
#         if(all(weight^exponent[l]==0)){
#           glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lambda=99e99,alpha=alpha.two)
#         } else {
#           glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lower.limits=rest.int$lower.limits,upper.limits=rest.int$upper.limits,penalty.factor=1/(weight^exponent[l]),alpha=alpha.two)
#         }
#         y_hat[[i]][[l]][cond,] <- predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=glm.two.ext[[i]][[l]]$lambda,type="response")
#       }
#     }
#   }
#   
#   metric <- list()
#   id.min <- lambda.min <- rep(x=NA,times=q)
#   exp.min <- rep(x=NA,times=q)
#   for(i in seq_len(q)){
#     metric[[i]] <- list()
#     for(l in seq_along(exponent)){
#       metric[[i]][[l]] <- apply(X=y_hat[[i]][[l]],MARGIN=2,FUN=function(x)
#         calc.metric(y=y[[i]],y_hat=x,family=family))
#     }
#     min <- sapply(metric[[i]],min)
#     tryCatch(expr=graphics::plot(x=exponent,y=min,type="o"),error=function(x) NULL)
#     id.exp <- which.min(min)
#     exp.min[i] <- exponent[id.exp]
#     id.min[i] <- which.min(metric[[i]][[id.exp]])
#     lambda.min[i] <- glm.two.ext[[i]][[id.exp]]$lambda[id.min[i]]
#     glm.two.ext[[i]] <- glm.two.ext[[i]][[id.exp]]
#   }
#   
#   list <- list(glm.one=glm.one.ext,glm.two=glm.two.ext,exp.min=exp.min,lambda.min=lambda.min,info=data.frame(p=p,q=q))
#   class(list) <- "sparselink"
#   return(list)
# }
# 
# 
# # This was a working version (2024-07-23).
# # (tune weighting of target data and source data)
# glm.both <- function(x,y,family,alpha.one=0.95,alpha.two=1){
#   warning("remove seed")
#   set.seed(1)
#   
#   n <- sapply(X=y,FUN=base::length)
#   p <- ncol(x[[1]])
#   q <- length(x)
#   
#   nfolds <- 10
#   foldid <- make.folds.trans(y=y,family=family,nfolds=nfolds)
#   
#   # overall fit
#   glm.one.ext <- glm.two.ext <- list()
#   alpha.ext <- rep(x=NA,times=q)
#   beta.ext <- pval.ext.inc <- pval.ext.dec <- matrix(data=NA,nrow=p,ncol=q)
#   for(i in seq_len(q)){
#     glm.one.ext[[i]] <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],family=family,alpha=alpha.one) # use fixed foldid?
#     temp <- coef(object=glm.one.ext[[i]],s="lambda.min")
#     alpha.ext[i] <- temp[1]
#     beta.ext[,i] <- temp[-1]
#   }
#   
#   cands <- seq(from=0,to=1,by=0.2) # 
#   weight <- expand.grid(source=cands,target=cands) # only for 2nd and 3rd schemes
#   
#   #weight <- data.frame(source=seq(from=0,to=1,by=0.1)) # only for 1st scheme
#   
#   for(i in seq_len(q)){
#     rest.ext <- comb_split(coef=beta.ext,weight=n,id=i)
#     glm.two.ext[[i]] <- list()
#     for(l in seq_len(nrow(weight))){
#       #pf <- 1/(weight$source[l]*rest.ext$weight.source + weight$target[l]*rest.ext$weight.target)
#       #pf <- 1/(weight$source[l]*rest.ext$weight.source + (1-weight$source[l])*rest.ext$weight.target)
#       #pf <- 1/(weight$source[l]*rest.ext$weight.source + weight$target[l]*rest.ext$weight.target) + (weight$source[l]==0 & weight$target[l]==0)
#       pf <- 1/(rest.ext$weight.source^weight$source[l] + rest.ext$weight.target^weight$target[l]) #- (weight$source[l]==0) - (weight$target[l]==0) + (weight$source[l]==0 & weight$target[l]==0))
#       if(all(pf==Inf)){
#         glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lambda=99e99,alpha=alpha.two)
#       } else {
#         glm.two.ext[[i]][[l]] <- glmnet::glmnet(x=cbind(x[[i]],x[[i]]),y=y[[i]],family=family,lower.limits=rest.ext$lower.limits,upper.limits=rest.ext$upper.limits,penalty.factor=pf,alpha=alpha.two)
#       }
#     }
#   }
#   
#   # cross-validation
#   y_hat <- list()
#   for(i in seq_len(q)){
#     y_hat[[i]] <- list()
#     for(l in seq_len(nrow(weight))){
#       nlambda <- length(glm.two.ext[[i]][[l]]$lambda)
#       y_hat[[i]][[l]] <- matrix(data=NA,nrow=n[i],ncol=nlambda)
#     }
#   }
#   
#   # Changed order of loops on 2024-07-22.
#   
#   #for(k in seq_len(nfolds)){ # original
#   #  alpha.int <- alpha.ext # original
#   #  beta.int <- beta.ext # original
#   #  for(i in seq_len(q)){ # original
#   
#   for(i in seq_len(q)){ # trial
#     alpha.int <- alpha.ext # trial
#     beta.int <- beta.ext # trial
#     for(k in seq_len(nfolds)){ # trial
#       
#       cond <- foldid[[i]]==k
#       glm.one.int <- glmnet::glmnet(x=x[[i]][!cond,],y=y[[i]][!cond],family=family,alpha=alpha.one)
#       temp <- coef(object=glm.one.int,s=glm.one.ext[[i]]$lambda.min)
#       alpha.int[i] <- temp[1]
#       beta.int[,i] <- temp[-1]
#       
#       rest.int <- comb_split(coef=beta.int,weight=n,id=i)
#       
#       for(l in seq_len(nrow(weight))){
#         #pf <- 1/(weight$source[l]*rest.int$weight.source + weight$target[l]*rest.int$weight.target)
#         #pf <- 1/(weight$source[l]*rest.int$weight.source + (1-weight$source[l])*rest.int$weight.target)
#         #pf <- 1/(weight$source[l]*rest.int$weight.source + weight$target[l]*rest.int$weight.target) + (weight$source[l]==0 & weight$target[l]==0)
#         pf <- 1/(rest.int$weight.source^weight$source[l] + rest.int$weight.target^weight$target[l]) #- (weight$source[l]==0) - (weight$target[l]==0) + (weight$source[l]==0 & weight$target[l]==0))
#         if(all(pf==Inf)){
#           glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lambda=99e99,alpha=alpha.two)
#         } else {
#           glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family,lower.limits=rest.int$lower.limits,upper.limits=rest.int$upper.limits,penalty.factor=pf,alpha=alpha.two)
#         }
#         y_hat[[i]][[l]][cond,] <- predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=glm.two.ext[[i]][[l]]$lambda,type="response")
#       }
#     }
#   }
#   
#   metric <- list()
#   id.min <- lambda.min <- rep(x=NA,times=q)
#   for(i in seq_len(q)){
#     metric[[i]] <- list()
#     for(l in seq_len(nrow(weight))){
#       metric[[i]][[l]] <- apply(X=y_hat[[i]][[l]],MARGIN=2,FUN=function(x)
#         calc.metric(y=y[[i]],y_hat=x,family=family))
#     }
#     min <- sapply(metric[[i]],min)
#     #tryCatch(expr=graphics::plot(x=seq_len(nrow(weight)),y=log(min),type="o"),error=function(x) NULL)
#     #tryCatch(expr=scatterplot3d::scatterplot3d(x=weight$source,y=weight$target,z=min,type="h",color="red",pch=16),error=function(x) NULL)
#     tryCatch(expr=graphics::plot(x=weight$source,y=log(min),type="o"),error=function(x) NULL)
#     
#     id.exp <- which.min(min)
#     id.min[i] <- which.min(metric[[i]][[id.exp]])
#     lambda.min[i] <- glm.two.ext[[i]][[id.exp]]$lambda[id.min[i]]
#     glm.two.ext[[i]] <- glm.two.ext[[i]][[id.exp]]
#   }
#   
#   list <- list(glm.one=glm.one.ext,glm.two=glm.two.ext,lambda.min=lambda.min,info=data.frame(p=p,q=q))
#   class(list) <- "sparselink"
#   return(list)
# }
# 
# # simulate toy data for multi-task and transfer learning
# if(FALSE){
#   
#   # multi-task learning  
#   n <- 100
#   p <- 200
#   q <- 3
#   x <- matrix(data=stats::rnorm(n*p),nrow=n,ncol=p)
#   beta.com <- stats::rbinom(n=p,size=1,prob=0.1)*stats::rnorm(n=p)
#   beta.sep <- matrix(data=NA,nrow=p,ncol=q)
#   for(i in seq_len(q)){
#     beta.sep[,i] <- beta.com + stats::rbinom(n=p,size=1,prob=0.05)*stats::rnorm(n=p)
#   }
#   y <- x %*% beta.sep + stats::rnorm(n=n*q)
#   family <- "gaussian"
#   
#   # transfer learning
#   n <- c(50,100,200)
#   p <- 200
#   q <- length(n)
#   x <- list()
#   beta.com <- stats::rbinom(n=p,size=1,prob=0.1)*stats::rnorm(n=p)
#   beta.sep <- list()
#   y <- list()
#   for(i in seq_len(q)){
#     x[[i]] <- matrix(data=stats::rnorm(n=n[i]*p),nrow=n[i],ncol=p)
#     beta.sep[[i]] <- beta.com + stats::rbinom(n=p,size=1,prob=0.05)*stats::rnorm(n=p)
#     y[[i]] <- x[[i]] %*% beta.sep[[i]] + stats::rnorm(n=n[i])
#   }
#   family <- "gaussian"
#   
#   # comparison
#   object <- sparselink(x=x,y=y,family="gaussian",type="exp")
#   
# }

# This is the current working version (2024-07-23).
# (multi-task learning and transfer learning)
# (using all folds for support problems in the case of transfer learning, using only training folds for all problems in the case of multi-target learning)


#' @title Construct penalty factors
#' 
#' @description
#' Uses internal and external weights (for negative and positive effects)
#' as well as internal and external exponents/factors for these weights
#' to construct penalty factors.
#' 
#' @param w_int internal weights:
#' numeric vector of length p with non-negative entries
#' @param w_ext external weights:
#' numeric vector of length p with non-negative entries
#' @param v_int exponent or factor for internal weights:
#' non-negative scalar
#' @param v_ext exponent or factor for external weights:
#' non-negative scalar
#' @param type character "geo", "exp", "rem" or "ari"
#' (with or without addition of ".con")
#' 
construct_pf <- function(w_int,w_ext,v_int,v_ext,type){
  if(type %in% c("geo","geo.con")){
    pf <- 1/((w_int^v_int)*(w_ext^v_ext))
  } else if(type %in% c("exp","exp.con")){
    pf <- 1/((w_int^v_int)+(w_ext^v_ext))
  } else if(type %in% c("rem","rem.con")){ 
    pf <- 1/((w_int^v_int)+(w_ext^v_ext)-1*(v_int==0)-1*(v_ext==0))
  } else if(type %in% c("ari","ari.con")){
    pf <- 1/(v_int*w_int + v_ext*w_ext)
  }
  if(any(pf<0)){stop("negative pf")}
  return(pf)
}

#' @title Visualise metric that depends on two parameters
#' @export
#' @keywords internal
#' 
#' @description 
#' Displays values in y in a grey scale (white=lowest, black=highest),
#' for different combinations of the two variables in x.
#' The lowest value is indicated by a red cross,
#' and the lowest value on the diagonal is indicated by a red circle.
#'
#' @param x list with slots source and target
#' @param y numeric vector
#' 
#' @examples
#' values <- seq(from=0,to=1,by=0.2)
#' x <- expand.grid(source=values,target=values)
#' y <- stats::rexp(n=length(values)*length(values))
#' plotWeight(x=x,y=y)
#' 
plotWeight <- function(x,y){
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

#' @title Sparse regression for related problems
#' @export
#' @param x n x p matrix (multi-task learning) 
#' or list of n_k x p matrices (transfer learning)
#' @param y n x q matrix (multi-task learning)
#' or list of n_k-dimensional vectors (transfer learning)
#' @param family character "gaussian" or "binomial"
#' @param alpha.init elastic net mixing parameter for initial regressions,
#' default: 0.95 (lasso-like elastic net)
#' @param alpha elastic net mixing parameter of final regressions,
#' default: 1 (lasso)
#' @param nfolds number of cross-validation folds
#' @param type character
#' 
#' @examples
#' # multi-task learning
#' n <- 100
#' p <- 50
#' q <- 3
#' family <- "gaussian"
#' x <- matrix(data=rnorm(n=n*p),nrow=n,ncol=p)
#' y <- matrix(data=rnorm(n*q),nrow=n,ncol=q)
#' object <- sparselink(x=x,y=y,family=family)
#' 
sparselink <- function(x,y,family,alpha.init=0.95,alpha=1,type="exp",nfolds=10){ # was alpha.init=0.95 and alpha=1
  
  alpha.one <- alpha.init
  alpha.two <- alpha
  
  #warning("remove seed")
  #set.seed(1)
  cat("type=",type,"\n")
  
  #trial <- TRUE # was FALSE
  #mode <- ""
  
  if(is.matrix(y) & is.matrix(x)){
    message("mode: multi-target learning")
    p <- ncol(x)
    q <- ncol(y)
    n <- rep(x=nrow(x),times=q)
    foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
    y <- apply(y,2,function(x) x,simplify=FALSE)
    x <- replicate(n=q,expr=x,simplify=FALSE)
    # Write function for creating fold identifiers (balancing unique rows of binary columns of target matrix).
    #foldid <- sample(rep(x=seq_len(nfolds),length.out=n[1]))
    foldid <- replicate(n=q,expr=foldid,simplify=FALSE)
    mode <- "multiple"
  } else if(is.list(y) & is.list(x)){
    message("mode: transfer learning")
    n <- sapply(X=y,FUN=base::length)
    p <- ncol(x[[1]])
    q <- length(x)
    foldid <- make.folds.trans(y=y,family=family,nfolds=nfolds)
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
  
  #warning("remove temporary line!")
  #mode <- "multiple"
  
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
  
  #weight <- c(seq(from=0,to=0.1,by=0.025),0.15,seq(from=0.2,to=0.8,by=0.1),0.85,seq(from=0.9,to=1.0,by=0.025))
  #warning("remove next line")
  #weight[weight==0] <- 0.01
  #weight[weight==1] <- 0.99
  #weight <- seq(from=0,to=0.5,by=0.1)
  
  if(type %in% c("geo","exp","rem","ari")){
    cands <- seq(from=0,to=1,by=0.2)
    weight <- expand.grid(source=cands,target=cands)
  } else if(type %in% c("geo.con","exp.con","rem.con","ari.con")){
    cands <- seq(from=0,to=1,by=0.05)
    weight <- data.frame(source=cands,target=1-cands)
  } else {
    stop(paste0("Invalid type=",type))
  }
  
  for(i in seq_len(q)){
    #rest.ext <- comb_split(coef=beta.ext,weight=n,id=i)
    rest.ext <- comb_split_trial(coef=beta.ext,id=i)
    glm.two.ext[[i]] <- list()
    for(l in seq_len(nrow(weight))){
      
      # if(trial){
      #   pf <- 1/(rest.ext$weight.source^weight$source[l]*rest.ext$weight.target^weight$target[l]) # trial!
      #   #pf <- 1/(rest.ext$weight.source^weight$source[l]+rest.ext$weight.target^weight$target[l]) #- 1*(weight$source[l]==0) - 1*(weight$target[l]==0)) # original #-(weight[l]==0|weight[l]==1)
      #   #pf <- 1/(rest.ext$weight.source^weight[l]+rest.ext$weight.target^(1-weight[l])) # -(weight[l]==0|weight[l]==1)
      # } else {
      #   pf <- 1/(weight$source[l]*rest.ext$weight.source + weight$target[l]*rest.ext$weight.target)
      # }
      # if(any(pf<0)){stop("negative pf")}
      
      pf <- construct_pf(w_int=rest.ext$weight.target,w_ext=rest.ext$weight.source,v_int=weight$target[l],v_ext=weight$source[l],type=type)
      
      
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
      
      #rest.int <- comb_split(coef=beta.int,weight=n,id=i)
      rest.int <- comb_split_trial(coef=beta.int,id=i)
      
      for(l in seq_len(nrow(weight))){
        
        # if(trial){
        #   pf <- 1/(rest.int$weight.source^weight$source[l]*rest.int$weight.target^weight$target[l]) # trial!
        #   #pf <- 1/(rest.int$weight.source^weight$source[l]+rest.int$weight.target^weight$target[l]) #- 1*(weight$source[l]==0) - 1*(weight$target[l]==0)) # original
        #   #pf <- 1/(rest.int$weight.source^weight[l]+rest.int$weight.target^(1-weight[l])) # -(weight[l]==0|weight[l]==1)
        # } else {
        #   pf <- 1/(weight$source[l]*rest.int$weight.source + weight$target[l]*rest.int$weight.target)
        # }
        # if(any(pf<0)){stop("negative pf")}
        
        pf <- construct_pf(w_int=rest.int$weight.target,w_ext=rest.int$weight.source,v_int=weight$target[l],v_ext=weight$source[l],type=type)
        
        if(all(pf==Inf)){
          glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family[i],lambda=99e99,alpha=alpha.two)
        } else {
          glm.two.int <- glmnet::glmnet(x=cbind(x[[i]],x[[i]])[!cond,],y=y[[i]][!cond],family=family[i],lower.limits=rest.int$lower.limits,upper.limits=rest.int$upper.limits,penalty.factor=pf,alpha=alpha.two)
        }
        y_hat[[i]][[l]][cond,] <- stats::predict(object=glm.two.int,newx=cbind(x[[i]],x[[i]])[cond,],s=glm.two.ext[[i]][[l]]$lambda,type="response")
      }
    }
  }
  
  # # original (until 2024-11-27)
  # metric <- list()
  # id.min <- lambda.min <- rep(x=NA,times=q)
  # for(i in seq_len(q)){
  #   metric[[i]] <- list()
  #   for(l in seq_len(nrow(weight))){
  #     metric[[i]][[l]] <- apply(X=y_hat[[i]][[l]],MARGIN=2,FUN=function(x)
  #       calc.metric(y=y[[i]],y_hat=x,family=family[i]))
  #   }
  #   min <- sapply(X=metric[[i]],FUN=min)
  #   tryCatch(expr=plotWeight(x=weight,y=log(min)),error=function(x) NULL)
  #   id.exp <- which.min(min)
  #   id.min[i] <- which.min(metric[[i]][[id.exp]])
  #   lambda.min[i] <- glm.two.ext[[i]][[id.exp]]$lambda[id.min[i]]
  #   glm.two.ext[[i]] <- glm.two.ext[[i]][[id.exp]]
  # }
  
  metric <- list()
  cvm.min <- lambda.ind <- lambda.min <- matrix(data=NA,nrow=nrow(weight),ncol=q)
  for(i in seq_len(q)){
    metric[[i]] <- list()
    for(l in seq_len(nrow(weight))){
      metric[[i]][[l]] <- apply(X=y_hat[[i]][[l]],MARGIN=2,FUN=function(x)
      calc.metric(y=y[[i]],y_hat=x,family=family[i]))
      lambda.ind[l,i] <- which.min(metric[[i]][[l]])
      lambda.min[l,i] <- glm.two.ext[[i]][[l]]$lambda[lambda.ind[l,i]]
    }
    cvm.min[,i] <- sapply(X=metric[[i]],FUN=min)
  }
  weight.ind <- apply(X=cvm.min,MARGIN=2,FUN=which.min)
  weight.min <- weight[weight.ind,]
  
  list <- list(glm.one=glm.one.ext,glm.two=glm.two.ext,weight=weight,weight.ind=weight.ind,weight.min=weight.min,lambda.min=lambda.min,info=data.frame(p=p,q=q,mode=mode,family=paste0(family,collapse=", ")))
  class(list) <- "sparselink"
  return(list)
}

# 
# if(FALSE){
#   p <- 50
#   x <- stats::rbinom(n=p,size=1,prob=0.4)*stats::runif(n=p)  
#   y <- stats::rbinom(n=p,size=1,prob=0.4)*stats::runif(n=p)
#   #wx <- 0.0
#   #wy <- 0.0
#   #z <- wx*x + wy * y
#   #z <- x^wx + y^wy - (wx==0) - (wy==0) + (wx==0 & wy==0)
#   #z <- z/sum(z)
#   #graphics::par(mfrow=c(1,2))
#   w <- 0.0
#   #z <- w*x + (1-w)*y
#   z <- x^w + y^(1-w) - (w==0|w==1)
#   ymax <- max(z)
#   graphics::plot(x=x,y=z,ylim=c(0,ymax))
#   graphics::plot(x=y,y=z,ylim=c(0,ymax))
#   
# }


if(FALSE){
  # Search for global variables!
  fun <- objects(all.names=TRUE,pattern="\\.")
  fun <- objects()
  for(i in seq_along(fun)){
    if(fun[i] %in% c("fun","path","i")){next}
    if(fun[i]==".Random.seed"){next}
    var <- codetools::findGlobals(fun=fun[i],merge=FALSE)$variables
    if(length(var)>0){
      warning(paste0("Global variable(s) in function ",fun[i],": ",paste0(var,collapse=", ")))
    }
  }
}
