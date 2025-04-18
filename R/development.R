
fuse <- function(x,y,mode){
  slot <- list(sep=list(),com=numeric())
  list <- list(y=slot,x=slot)
  #--- standardise data (MTL and TL) ---
  for(i in seq_along(x)){
    list$x$sep[[i]] <- scale(x[[i]])
  }
  for(i in seq_along(y)){
    list$y$sep[[i]] <- scale(y[[i]])
  }
  #--- average (MTL) or concatenate (TL) data ---
  if(mode=="multiple"){
    list$y$com <- rowSums(sapply(X=y,FUN=scale))
    list$x$com <- x[[1]]
  } else if(mode=="transfer"){
    list$y$com <- unlist(lapply(X=y,FUN=scale))
    list$x$com <- do.call(what="rbind",args=x)
  }
  return(list)
}

# alpha=NA returns Spearman's correlation coefficients
init.coef <- function(list,alpha=0.95,lambda.sep=NULL,lambda.com=NULL,trial=FALSE){
  if(is.null(lambda.sep)!=is.null(lambda.com)){stop()}
  cv <- is.null(lambda.sep) & is.null(lambda.com)
  q <- length(list$x$sep)
  glm.sep <- coef.sep <- list()
  if(!is.na(alpha) & cv){
    lambda.sep <- numeric() 
  }
  #--- group lasso (currently only for Gaussian MTL)
  if(trial){
    cat("trial","\n")
    y_temp <- sapply(list$y$sep,function(x) x)
    x_temp <- list$x$sep[[1]]
    if(cv){
      object <- glmnet::cv.glmnet(x=x_temp,y=y_temp,family="mgaussian",alpha=alpha)
      coef <- stats::coef(object,s="lambda.min")
      lambda.com <- lambda.sep <- object$lambda.min
    } else {
      object <- glmnet::glmnet(x=x_temp,y=y_temp,family="mgaussian",alpha=alpha)
      coef <- stats::coef(object=object,s=lambda.com)
    }
    coef.sep <- lapply(coef,function(x) x[-1])
    coef.com <- do.call(what="+",args=coef.sep)
  }
  
  if(!trial){
  #--- separate models ---
  for(i in seq_len(q)){
    if(is.na(alpha)){
      coef.sep[[i]] <- as.numeric(stats::cor(x=list$x$sep[[i]],y=list$y$sep[[i]],method="spearman"))
    } else if(cv){
      glm.sep[[i]] <- glmnet::cv.glmnet(x=list$x$sep[[i]],y=list$y$sep[[i]],alpha=alpha)
      coef.sep[[i]] <- stats::coef(glm.sep[[i]],s="lambda.min")[-1]
      lambda.sep[i] <- glm.sep[[i]]$lambda.min
    } else {
      glm.sep[[i]] <- glmnet::glmnet(x=list$x$sep[[i]],y=list$y$sep[[i]],alpha=alpha)
      coef.sep[[i]] <- stats::coef(glm.sep[[i]],s=lambda.sep[i])[-1]
    }
  }
  #--- common model ---
  if(is.na(alpha)){
    n <- sapply(X=list$x$sep,FUN=nrow)
    coef.com <- sapply(X=coef.sep,function(x) x) %*% n/sum(n) #
    # Use Fisher transform instead!
    warning("Use Fisher transform.") # no 
  } else if(cv){
    glm.com <- glmnet::cv.glmnet(x=list$x$com,y=list$y$com,alpha=alpha)
    coef.com <- stats::coef(glm.com,s="lambda.min")[-1]
    lambda.com <- glm.com$lambda.min
  } else {
    glm.com <- glmnet::glmnet(x=list$x$com,y=list$y$com,alpha=alpha)
    coef.com <- stats::coef(glm.com,s=lambda.com)[-1]
  }
  }  
    
  list <- list(com=coef.com,sep=coef.sep,lambda.sep=lambda.sep,lambda.com=lambda.com)
  return(list)
}

penfac <- function(sep,com,exp.sep,exp.com){ # prop, 
  # combining with proportion
  #positive <- (1-prop)*pmax(0,sep)+prop*pmax(0,com)
  #negative <- -(1-prop)*pmin(0,sep)-prop*pmin(0,com)
  # combining with exponents
  positive <- pmax(0,sep)^exp.sep + pmax(0,com)^exp.com
  negative <- (-pmin(0,sep))^exp.sep + (-pmin(0,com))^exp.com
  pf <- 1/c(positive,negative)
  pf[pf==-Inf] <- Inf
  return(pf)
}

devel_old <- function(x,y,family="gaussian",alpha.init=0.95,alpha=1,nfolds=10,trial=FALSE){
  if(any(family!="gaussian")){stop("not implemented")}
  #if(alpha!=1){stop("not implemented")}
  
  if(is.matrix(y) & is.matrix(x)){
    message("mode: multi-target learning")
    p <- ncol(x)
    q <- ncol(y)
    n <- rep(x=nrow(x),times=q)
    foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
    y <- apply(X=y,MARGIN=2,FUN=function(x) x,simplify=FALSE)
    x <- replicate(n=q,expr=x,simplify=FALSE)
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
    family <- rep(x=family,times=q)
  }
  
  list.ext <- fuse(x=x,y=y,mode=mode)
  init.ext <- init.coef(list=list.ext,alpha=alpha.init,trial=trial)

  #ncand <- 11
  #prop <- seq(from=0,to=1,length.out=ncand)
  exp <- c(0.0,0.2,0.5,1.0,2.0,5.0) # flexible
  #exp <- seq(from=0,to=1,by=0.2) # unit interval
  grid <- expand.grid(sep=exp,com=exp)
  ncand <- nrow(grid)
  model.ext <- list()
  for(j in seq_len(q)){
    model.ext[[j]] <- list()
    for(k in seq_len(ncand)){
      pf <- penfac(sep=init.ext$sep[[j]],com=init.ext$com,exp.sep=grid$sep[k],exp.com=grid$com[k]) # prop=grid$prop[k]
      if(all(is.infinite(pf))){
        model.ext[[j]][[k]] <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]]),y=y[[j]],family=family[j],lambda=99e99,alpha=alpha)
      } else {
        model.ext[[j]][[k]] <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]]),y=y[[j]],family=family[j],lower.limits=0,penalty.factor=pf,alpha=alpha)
      }
    }
  }
  
  #--- initialise matrices ---
  pred <- list()
  for(j in seq_len(q)){
    pred[[j]] <- list()
    for(k in seq_len(ncand)){
      pred[[j]][[k]] <- matrix(data=NA,nrow=n[j],ncol=length(model.ext[[j]][[k]]$lambda))
    }
  }
  
  #--- cross-validation ---
  for(i in seq_len(nfolds)){
    x.train <- y.train <- x.test <- list()
    for(j in seq_len(q)){
      cond <- foldid[[j]]==i
      x.train[[j]] <- x[[j]][!cond,,drop=FALSE]
      y.train[[j]] <- y[[j]][!cond]
      x.test[[j]] <- x[[j]][cond,,drop=FALSE]
    }
    list.int <- fuse(x=x.train,y=y.train,mode=mode)
    init.int <- init.coef(list=list.int,lambda.sep=init.ext$lambda.sep,lambda.com=init.ext$lambda.com,alpha=alpha.init,trial=trial)

    for(j in seq_len(q)){
      cond <- foldid[[j]]==i
      for(k in seq_len(ncand)){
        pf <- penfac(sep=init.int$sep[[j]],com=init.int$com,exp.sep=grid$sep[k],exp.com=grid$com[k]) # prop=grid$prop[k],
        if(all(is.infinite(pf))){
          model.int <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]])[!cond,],y=y[[j]][!cond],family=family[j],lambda=99e99,alpha=alpha)
        } else {
          model.int <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]])[!cond,],y=y[[j]][!cond],family=family[j],lower.limits=0,penalty.factor=pf,alpha=alpha)
        }
        pred[[j]][[k]][cond,] <- stats::predict(object=model.int,newx=cbind(x[[j]],-x[[j]])[cond,],s=model.ext[[j]][[k]]$lambda)
      }
    }
  }
  
  #--- calculate MSE ---
  mse <- list()
  id.grid <- lambda.min <- numeric()
  for(j in seq_len(q)){
    mse[[j]] <- list()
    for(k in seq_len(ncand)){
      mse[[j]][[k]] <- apply(X=pred[[j]][[k]],MARGIN=2,FUN=function(x) mean((x-y[[j]])^2))
    }
    id.grid[j] <- which.min(sapply(mse[[j]],min))
    lambda.min[j] <- model.ext[[j]][[id.grid[j]]]$lambda[which.min(mse[[j]][[id.grid[j]]])]
    #graphics::plot(sapply(X=mse[[j]],FUN=min),type="o")
    temp <- grid
    names(temp) <- c("source","target")
    tryCatch(expr=plotWeight(x=temp,y=sapply(X=mse[[j]],FUN=min)),error=function(x) NULL)
  }
  
  list <- list(model=model.ext,id.grid=id.grid,lambda.min=lambda.min)
  class(list) <- "devel"
  return(list)
}

predict.devel <- function(object,newx){
  q <- length(object$model)
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- stats::predict(object=object$model[[i]][[object$id.grid[i]]],newx=cbind(newx,-newx),s=object$lambda.min[i])
  }
  return(y_hat)
}

coef.devel <- function(object){
  return(list(alpha=NA,beta=NA))
}


#---- correlation-based re-implementation ---

#cordev.init <- function(x,y,alpha.init=NA){
#  p <- ncol(x[[1]])
#  q <- length(x)
#  cor <- matrix(data=NA,nrow=p,ncol=q)
#  for(i in seq_len(q)){
#    cor[,i] <- stats::cor(y=y[[i]],x=x[[i]],method="spearman")
#  }
#  return(cor)
#}

cordev <- function(x,y,family="gaussian",alpha.init=NA,nfolds=10){
  
  if(is.matrix(y) & is.matrix(x)){
    message("mode: multi-target learning")
    p <- ncol(x)
    q <- ncol(y)
    n <- rep(x=nrow(x),times=q)
    foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
    y <- apply(X=y,MARGIN=2,FUN=function(x) x,simplify=FALSE)
    x <- replicate(n=q,expr=x,simplify=FALSE)
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
  
  cor.ext <- matrix(data=NA,nrow=p,ncol=q)
  for(i in seq_len(q)){
    cor.ext[,i] <- stats::cor(y=y[[i]],x=x[[i]],method="spearman")
  }
  #cor.ext <- cordev.init(x=x,y=y,alpha.init=alpha.init)
  #rel.ext <- stats::cor(cor.ext,method="spearman")
  
  cand <- c(0,0.2,0.5,1,2,5)
  #grid <- expand.grid(sep=cand,com=cand)
  grid <- expand.grid(com=seq(from=0,to=10,length.out=21))
  
  weight.ext <- list()
  #weight.ext$ind <- rbind(pmax(cor.ext,0),-pmin(cor.ext,0))
  # switch for Fisher-transform (next line, matters for TL; at least account for n)
  weight.ext$com <- c(rowMeans(pmax(cor.ext,0)),rowMeans(-pmin(cor.ext,0)))
  object.ext <- list()
  for(i in seq_len(q)){
    object.ext[[i]] <- list()
    for(j in seq_len(nrow(grid))){
      #weight <- rel.ext[i,]
      #weight[i] <- 0
      # Remove multiplication with %*% weight outside of pmax/pmin!!!
      #temp <- rbind(pmax(cor.ext %*% weight,0),-pmin(cor.ext %*% weight,0))
      #temp <- c(rowSums(pmax(cor.ext,0)),rowSums(-pmin(cor.ext,0)))
      #temp <- c(rowSums(pmax(t(t(cor.ext)*weight),0)),rowSums(-pmin(t(t(cor.ext)*weight),0)))
      #pf.ext <- 1/(weight.ext$ind[,i]^grid$sep[j]+temp^grid$com[j])
      pf.ext <- 1/(weight.ext$com^grid$com[j])
      object.ext[[i]][[j]] <- glmnet::glmnet(x=cbind(x[[i]],-x[[i]]),y=y[[i]],family=family,lower.limits=0,penalty.factor=pf.ext)
    }
  }
  
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- list()
    for(j in seq_len(nrow(grid))){
      y_hat[[i]][[j]] <- matrix(data=NA,nrow=n[i],ncol=length(object.ext[[i]][[j]]$lambda))
    }
  }
  
  # cross-validation
  for(k in seq_len(nfolds)){
    cond <- foldid[[i]]==k
    cor.int <- matrix(data=NA,nrow=p,ncol=q)
    for(i in seq_len(q)){
      cor.int[,i] <- stats::cor(y=y[[i]][!cond],x=x[[i]][!cond,],method="spearman")
    }
    #cor.int <- cordev.init(x=x[[i]][!cond,],y=y[[i]][!cond],alpha.init=alpha.init)
    #rel.int <- stats::cor(cor.int,method="spearman")
    weight.int <- list()
    #weight.int$ind <- rbind(pmax(cor.int,0),-pmin(cor.int,0))
    weight.int$com <- c(rowMeans(pmax(cor.int,0)),rowMeans(-pmin(cor.int,0)))
    object.int <- list()
    for(i in seq_len(q)){
      for(j in seq_len(nrow(grid))){
        #weight <- rel.int[i,]
        #weight[i] <- 0
        #temp <- rbind(pmax(cor.int %*% weight,0),-pmin(cor.int %*% weight,0))
        #temp <- c(rowSums(pmax(cor.int,0)),rowSums(-pmin(cor.int,0)))
        #temp <- c(rowSums(pmax(t(t(cor.int)*weight),0)),rowSums(-pmin(t(t(cor.int)*weight),0)))
        #pf.int <- 1/(weight.int$ind[,i]^grid$sep[j]+temp^grid$com[j])
        pf.int <- 1/(weight.int$com^grid$com[j])
        object.int <- glmnet::glmnet(x=cbind(x[[i]],-x[[i]])[!cond,],y=y[[i]][!cond],family=family,lower.limits=0,penalty.factor=pf.int)
        y_hat[[i]][[j]][cond] <- predict(object=object.int,newx=cbind(x[[i]],-x[[i]])[cond,],s=object.ext[[i]][[j]]$lambda,type="response") 
      }
    }
  }
  
  # metric
  metric <- list()
  lambda.min <- numeric()
  for(i in seq_len(q)){
    metric[[i]] <-  list()
    for(j in seq_len(nrow(grid))){
      metric[[i]][[j]] <- apply(y_hat[[i]][[j]],2,function(x) mean((y[[i]]-x)^2))
      #lambda.min[[i]][j] <- object.ext[[i]][[j]]$lambda[which.min(metric[[i]][[j]])]
    }
    id <- which.min(sapply(metric[[i]],min))
    object.ext[[i]] <- object.ext[[i]][[id]]
    metric[[i]] <- metric[[i]][[id]]
    lambda.min[i] <- object.ext[[i]]$lambda[which.min(metric[[i]])]
  }
  
  #grid.min <- lapply(metric,function(x) sapply(x,which.min))
  
  list <- list(model=object.ext,lambda.min=lambda.min)
  class(list) <- "cordev"
  return(list)
}

predict.cordev <- function(object,newx){
  q <- length(object$model)
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- stats::predict(object=object$model[[i]],newx=cbind(newx,-newx),s=object$lambda.min[i])
  }
  return(y_hat)
}

coef.cordev <- function(object){
  return(list(alpha=NA,beta=NA))
}

#--- group-lasso based re-implementation ---

group.init <- function(x,y,alpha.init,lambda.min=NULL){
  if(is.null(lambda.min)){
    object <- glmnet::cv.glmnet(x=x,y=y,family="mgaussian",alpha.init=alpha.init)
    coef <- sapply(stats::coef(object=object,s="lambda.min"),function(x) x[-1])
    lambda.min <- object$lambda.min
  } else {
    object <- glmnet::glmnet(x=x,y=y,family="mgaussian",alpha.init=alpha.init)
    coef <- sapply(stats::coef(object=object,s=lambda.min),function(x) x[-1])
  }
  sep <- rbind(pmax(coef,0),-pmin(coef,0))
  com <- c(rowSums(pmax(coef,0)),rowSums(-pmin(coef,0)))
  list <- list(sep=sep,com=com,lambda.min=lambda.min)
  return(list)
}

group.devel <- function(x,y,family="gaussian",nfolds=10,alpha=1,alpha.init=0.95){
  if(alpha!=1){warning("arg not implemented")}
  
  p <- ncol(x)
  q <- ncol(y)
  n <- rep(x=nrow(x),times=q)
  foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)

  init.ext <- group.init(x=x,y=y,alpha.init=alpha.init)
  
  #graphics::plot(x=beta,y=init.ext$com[1:p]-init.ext$com[(p+1):(2*p)])

  #cand <- c(0,0.2,0.5,1,2,5)
  #grid <- expand.grid(sep=cand,com=cand)
  cand <- seq(from=0,to=5,length.out=11)
  grid <- data.frame(w=cand)
  
  object.ext <- list()
  for(i in seq_len(q)){
    object.ext[[i]] <- list()
    for(j in seq_len(nrow(grid))){
      #pf.ext <- rep(1,times=2*p)# remove this line
      #pf.ext <- 1/(init.ext$sep[,i]^grid$sep[j]+init.ext$com^grid$com[j])
      #pf.ext <- 1/(init.ext$sep[,i]*(1-grid$w[j])+init.ext$com*grid$w[j])
      #pf.ext <- 1/(init.ext$com^grid$w[j])
      pf.ext <- 1/(init.ext$sep[,i]^grid$w[j])
      object.ext[[i]][[j]] <- glmnet::glmnet(x=cbind(x,-x),y=y[,i],penalty.factor=pf.ext,lower.limits=0)
      #beta <- coef(object.ext[[i]][[j]],s=0.1)[-1]
      #plot(x=1/pf.ext,y=beta)
    }; rm(j)
  }; rm(i)
  
  y_hat <- list()
  for(i in seq_len(q)){
    y_hat[[i]] <- list()
    for(j in seq_len(nrow(grid))){
      y_hat[[i]][[j]] <- matrix(data=NA,nrow=n,ncol=length(object.ext[[i]][[j]]$lambda))
    }; rm(j)
  }; rm(i)
  
  for(k in seq_len(nfolds)){
    init.int <- group.init(x=x[foldid!=k,],y=y[foldid!=k,],alpha.init=alpha.init,lambda.min=init.ext$lambda.min)
    for(i in seq_len(q)){
      for(j in seq_len(nrow(grid))){
        #pf.int <- rep(x=1,times=2*p) # remove this line
        #pf.int <- 1/(init.int$sep[,i]^grid$sep[j]+init.int$com^grid$com[j])
        #pf.int <- 1/(init.int$sep[,i]*(1-grid$w[j])+init.int$com*grid$w[j])
        #pf.int <- 1/(init.int$com^grid$w[j])
        pf.int <- 1/(init.int$sep[,i]^grid$w[j])
        object.int <- glmnet::glmnet(x=cbind(x,-x)[foldid!=k,],y=y[foldid!=k,i],penalty.factor=pf.int,lower.limits=0)
        y_hat[[i]][[j]][foldid==k,] <- stats::predict(object=object.int,newx=cbind(x,-x)[foldid==k,],s=object.ext[[i]][[j]]$lambda,type="response")
      }
    }
  }
  
  mse <- list()
  for(i in seq_len(q)){
    mse[[i]] <- list()
    for(j in seq_len(nrow(grid))){
      mse[[i]][[j]] <- apply(X=y_hat[[i]][[j]],MARGIN=2,FUN=function(x) mean((x-y[,i])^2))
    }
  }; rm(i)
  
  
  model <- list()
  lambda.min <- numeric()
  for(i in seq_len(q)){
    mse.min <- sapply(X=mse[[i]],FUN=min)
    id.grid <- which.min(mse.min)
    lambda.min[i] <- object.ext[[i]][[id.grid]]$lambda[which.min(mse[[i]][[id.grid]])]
    model[[i]] <- object.ext[[i]][[id.grid]]
  }; rm(i)
  
  list <- list(model=model,lambda.min=lambda.min)
  class(list) <- "cordev"
  return(list)
}


#- - - - - - - - - - - - - - - - - - -
#---- simplified re-implementation ----
#- - - - - - - - - - - - - - - - - - -

# alpha=NA returns Spearman's correlation coefficients
init.coef <- function(x,y,alpha=0.95,lambda=NULL){
  cv <- is.null(lambda)
  q <- length(x)
  glm <- coef <- list()
  if(!is.na(alpha) & cv){
    lambda <- numeric() 
  }
  #--- separate models ---
  for(i in seq_len(q)){
    if(is.na(alpha)){
      coef[[i]] <- as.numeric(stats::cor(x=x[[i]],y=y[[i]],method="spearman"))
    } else if(cv){
      glm[[i]] <- glmnet::cv.glmnet(x=x[[i]],y=y[[i]],alpha=alpha)
      coef[[i]] <- stats::coef(glm[[i]],s="lambda.min")[-1]
      lambda[i] <- glm[[i]]$lambda.min
    } else {
      glm[[i]] <- glmnet::glmnet(x=x[[i]],y=y[[i]],alpha=alpha)
      coef[[i]] <- stats::coef(glm[[i]],s=lambda[i])[-1]
    }
  }
  list <- list(coef=coef,lambda=lambda)
  return(list)
}

penfac <- function(coef,exp){
  coef <- sapply(coef,function(x) x)
  positive <- rowSums(pmax(coef,0))^exp
  negative <- rowSums(-pmin(coef,0))^exp
  pf <- 1/c(positive,negative)
  pf[pf==-Inf] <- Inf
  return(pf)
}

devel <- function(x,y,family="gaussian",alpha.init=0.95,alpha=1,nfolds=10){
  #if(any(family!="gaussian")){stop("not implemented")}
  #if(alpha!=1){stop("not implemented")}
  
  if(is.matrix(y) & is.matrix(x)){
    message("mode: multi-target learning")
    p <- ncol(x)
    q <- ncol(y)
    n <- rep(x=nrow(x),times=q)
    foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
    y <- apply(X=y,MARGIN=2,FUN=function(x) x,simplify=FALSE)
    x <- replicate(n=q,expr=x,simplify=FALSE)
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
    family <- rep(x=family,times=q)
  }
  
  init.ext <- init.coef(x=x,y=y,alpha=alpha.init)
  

  exp <- c(0.0,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,5.0) # flexible
  grid <- expand.grid(exp=exp)
  ncand <- nrow(grid)
  model.ext <- list()
  for(j in seq_len(q)){
    model.ext[[j]] <- list()
    for(k in seq_len(ncand)){
      pf <- penfac(coef=init.ext$coef,exp=grid$exp[k])
      if(all(is.infinite(pf))){
        model.ext[[j]][[k]] <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]]),y=y[[j]],family=family[j],lambda=99e99,alpha=alpha)
      } else {
        model.ext[[j]][[k]] <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]]),y=y[[j]],family=family[j],lower.limits=0,penalty.factor=pf,alpha=alpha)
      }
    }
  }
  
  #--- initialise matrices ---
  pred <- list()
  for(j in seq_len(q)){
    pred[[j]] <- list()
    for(k in seq_len(ncand)){
      pred[[j]][[k]] <- matrix(data=NA,nrow=n[j],ncol=length(model.ext[[j]][[k]]$lambda))
    }
  }
  
  #--- cross-validation ---
  for(i in seq_len(nfolds)){
    x.train <- y.train <- x.test <- list()
    for(j in seq_len(q)){
      cond <- foldid[[j]]==i
      x.train[[j]] <- x[[j]][!cond,,drop=FALSE]
      y.train[[j]] <- y[[j]][!cond]
      x.test[[j]] <- x[[j]][cond,,drop=FALSE]
    }
    init.int <- init.coef(x=x.train,y=y.train,lambda=init.ext$lambda,alpha=alpha.init)
    
    for(j in seq_len(q)){
      cond <- foldid[[j]]==i
      for(k in seq_len(ncand)){
        pf <- penfac(coef=init.int$coef,exp=grid$exp[k])
        if(all(is.infinite(pf))){
          model.int <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]])[!cond,],y=y[[j]][!cond],family=family[j],lambda=99e99,alpha=alpha)
        } else {
          model.int <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]])[!cond,],y=y[[j]][!cond],family=family[j],lower.limits=0,penalty.factor=pf,alpha=alpha)
        }
        pred[[j]][[k]][cond,] <- stats::predict(object=model.int,newx=cbind(x[[j]],-x[[j]])[cond,],s=model.ext[[j]][[k]]$lambda,type="response")
      }
    }
  }
  
  #--- calculate MSE ---
  mse <- list()
  id.grid <- lambda.min <- numeric()
  tryCatch(graphics::par(mfrow=c(1,q)))
  for(j in seq_len(q)){
    mse[[j]] <- list()
    for(k in seq_len(ncand)){
      #mse[[j]][[k]] <- apply(X=pred[[j]][[k]],MARGIN=2,FUN=function(x) mean((x-y[[j]])^2))
      mse[[j]][[k]] <- apply(X=pred[[j]][[k]],MARGIN=2,FUN=function(x) calc.metric(y=y[[j]],y_hat=x,family=family[j]))
    }
    id.grid[j] <- which.min(sapply(mse[[j]],min))
    lambda.min[j] <- model.ext[[j]][[id.grid[j]]]$lambda[which.min(mse[[j]][[id.grid[j]]])]
    tryCatch(expr=graphics::plot(x=grid$exp,y=sapply(X=mse[[j]],FUN=min),type="o"),error=function(x) NULL)
  }
  
  list <- list(model=model.ext,id.grid=id.grid,lambda.min=lambda.min)
  class(list) <- "devel"
  return(list)
}

predict.devel <- function(object,newx){
  q <- length(object$model)
  y_hat <- list()
  if(is.list(newx)){
    for(i in seq_len(q)){
      y_hat[[i]] <- stats::predict(object=object$model[[i]][[object$id.grid[i]]],newx=cbind(newx[[i]],-newx[[i]]),s=object$lambda.min[i],type="response")
    }
  } else {
    for(i in seq_len(q)){
      y_hat[[i]] <- stats::predict(object=object$model[[i]][[object$id.grid[i]]],newx=cbind(newx,-newx),s=object$lambda.min[i],type="response")
    }
  }
  return(y_hat)
}

coef.devel <- function(object){
  return(list(alpha=NA,beta=NA))
}

#- - - - - - - - - - - - - - - - - - 
#---- common support ----
#- - - - - - - - - - - - - - - - - - 

# CONTINUE HERE: Re-consider common support approach without group lasso.
# Add x on the left multiple times (to allow for common effects of x on y for all problems).
# Then cross-validate weight of common effects and of separate effects.

expand <- function(x,y=NULL,q=NULL){
  if(is.null(y)){
    yy <- NULL
  } else if(is.list(y)){
    yy <- do.call(what="c",args=y)
  } else {
    yy <- unlist(y)
  }
  if(is.list(x)){
    n <- sapply(x,nrow)
    p <- ncol(x[[1]])
    q <- length(x)
    intercept <- matrix(data=0,nrow=sum(n),ncol=q)
    xx <- matrix(data=0,nrow=sum(n),ncol=q*p)
    cumsum <- c(0,cumsum(n))
    for(i in seq_len(q)){
      rows <- (cumsum[i]+1):(cumsum[i+1])
      cols <- ((i-1)*p+1):(i*p)
      xx[rows,cols] <- x[[i]]
      intercept[rows,i] <- 1
    }
  } else {
    n <- nrow(x)
    p <- ncol(x)
    intercept <- matrix(data=0,nrow=q*n,ncol=q)
    xx <- matrix(data=0,nrow=q*n,ncol=q*p)
    for(i in seq_len(q)){
      rows <- ((i-1)*n+1):(i*n)
      cols <- ((i-1)*p+1):(i*p)
      xx[rows,cols] <- x
      intercept[rows,i] <- 1
    }
  }
  list <- list(xx=xx,yy=yy,intercept=intercept)
  return(list)
}


common_support <- function(x,y,family){
  if(is.matrix(y) & is.matrix(x)){
    message("mode: multi-target learning")
    p <- ncol(x)
    q <- ncol(y)
    n <- rep(x=nrow(x),times=q)
    foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
    y <- apply(X=y,MARGIN=2,FUN=function(x) x,simplify=FALSE)
    x <- replicate(n=q,expr=x,simplify=FALSE)
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
  #intercept <- rep(x=paste0("P",seq_len(q)),times=n)
  feature <- rep(x=seq_len(p),times=q)
  group <- c(feature,rep(p+1,times=q))
  pf_group <- rep(x=c(1,0),times=c(p,1))
  list <- expand(x=x,y=y)
  pf_sparse <- rep(x=c(0,1),times=c(p*q,q)) # to use standard group lasso (not sparse)
  object <- sparsegl::cv.sparsegl(x=cbind(list$xx,list$intercept),y=list$yy,group=group,pf_group=pf_group,pf_sparse=pf_sparse,family=family,intercept=FALSE)
  list <- list(object=object,q=q)
  class(list) <- "comsup"
  return(list)
}

predict.comsup <- function(object,newx){
  if(is.list(newx)){
    intercept <- rep(x=seq_along(x),times=sapply(x,nrow))
    list <- expand(x=newx)
  } else {
    intercept <- rep(x=seq_len(object$q),each=nrow(newx))
    list <- expand(x=newx,q=object$q)
  }
  predict(object=object$object,newx=cbind(list$xx,list$intercept),type="response",s="lambda.min")
}

#- - - - - - - - - - - - - - - - - - -
#---- exploratory simulation ----
#- - - - - - - - - - - - - - - - - - -

if(FALSE){
  #alpha.init <- 0
  metric <- list()
  for(k in 1:10){
    # simulate data
    n0 <- 100
    n1 <- 10000
    #n1 <- 0 # remove this line
    n <- n0 + n1
    p <- 200
    x <- matrix(data=stats::rnorm(n*p),nrow=n,ncol=p)
    beta <- stats::rbinom(n=p,size=1,prob=0.1)*stats::rnorm(n=p)
    eta <- as.numeric(x %*% beta)
    y1 <- 1+eta + 0.2*stats::rnorm(n=n,sd=sd(eta))
    y2 <- 2+eta + 0.5*stats::rnorm(n=n,sd=sd(eta))
    y3 <- 3+eta + 0.7*stats::rnorm(n=n,sd=sd(eta))
    #y1 <- eta + 0.2*stats::rnorm(n=n,sd=sd(eta))
    #y2 <- eta + 0.4*stats::rnorm(n=n,sd=sd(eta))
    #y3 <- eta + 0.6*stats::rnorm(n=n,sd=sd(eta))
    y <- cbind(y1,y2,y3)
    y <- scale(y)
    #y <- cbind(y1,y2)
    q <- ncol(y)
    fold <- rep(x=c(0,1),times=c(n0,n1))
    y_hat <- list()
    #--- intercept-only model ---
    y_hat$empty <- matrix(colMeans(y[fold==0,]),nrow=n1,ncol=q,byrow=TRUE)
    #--- standard lasso ---
    y_hat$lasso <- matrix(data=NA,nrow=n1,ncol=q)
    for(i in seq_len(q)){
      object <- glmnet::cv.glmnet(x=x[fold==0,],y=y[fold==0,i])
      y_hat$lasso[,i] <- predict(object=object,newx=x[fold==1,],s="lambda.min")
    }
    #---- mgaussian ---
    object <- glmnet::cv.glmnet(x=x[fold==0,],y=y[fold==0,],family="mgaussian")
    y_hat$mgaussian <- predict(object=object,newx=x[fold==1,],s="lambda.min")[,,1]
    #--- sparselink ---
    #warning("reset to alpha.init=0.95")
    object <- sparselink(x=x[fold==0,],y=y[fold==0,],family="gaussian",alpha.init=0.95)
    temp <- predict(object=object,newx=x[fold==1,])
    y_hat$sparselink <- do.call(what="cbind",args=temp)
    
    ##--- group lasso start ---
    #yy <- as.numeric(y)
    #xx <- rbind(x,x)
    #zz <- rep(c(0,1),each=n)
    #ff <- c(fold,fold)
    #xx_int <- cbind(xx,xx)
    #group <- c(1,rep(x=seq(from=2,to=p+1),times=2))
    ## define foldid! (putting all entries from the same sample in the same group)
    #test <- gglasso::cv.gglasso(x=cbind(zz,xx_int)[ff==0,],y=yy[ff==0],group=group,pf=c(0,rep(1,times=p)))
    ## CONTINUE HERE
    #temp <- predict(test,newx=cbind(zz,xx_int)[ff==1,])
    #y_hat$group <- matrix(temp,ncol=2)
    #--- group lasso end ---
    
    #--- correlation-based re-implementation of sparselink ---
    object <- cordev(x=x[fold==0,],y=y[fold==0,],family="gaussian")
    temp <- predict(object=object,newx=x[fold==1,])
    y_hat$cordev <- do.call(what="cbind",args=temp)
    
    #--- group-lasso based approach ---
    object <- group.devel(x=x[fold==0,],y=y[fold==0,],family="gaussian",alpha.init=0.95)
    temp <- predict(object=object,newx=x[fold==1,])
    y_hat$group <- do.call(what="cbind",args=temp)
    
    #--- development  ---
    object <- devel(x=x[fold==0,],y=y[fold==0,],alpha.init=0.95,family="gaussian")
    temp <- predict(object=object,newx=x[fold==1,])
    y_hat$devel <- do.call(what="cbind",args=temp)
    
    #--- common support ---
    object <- common_support(x=x[fold==0,],y=y[fold==0,],family="gaussian")
    temp <- predict(object=object,newx=x[fold==1,])
    y_hat$comsup <- matrix(data=temp,ncol=q)
    #image(matrix(coef(object$object)[1:(p*q)],nrow=p,ncol=q))
    
    #--- prediction error ---
    mse <- matrix(data=NA,nrow=length(y_hat),ncol=q,dimnames=list(names(y_hat),NULL))
    for(i in seq_along(y_hat)){
      for(j in seq_len(q)){
        mse[i,j] <- mean((y[fold==1,j]-y_hat[[i]][,j])^2)
      }
    }
    metric[[k]] <- mse
  }
  rowMeans(sapply(metric,function(x) colMeans(t(x)/x["empty",])))
  #rowMeans(do.call(what="cbind",args=metric))
  #Reduce(f="+",x=metric)
}

#object <- devel(x=x,y=y,family="gaussian")
#y_hat <- predict(object=object,newx=x)

