
fuse <- function(x,y,mode){
  slot <- list(sep=list(),com=numeric())
  list <- list(y=slot,x=slot)
  # standardise data (MTL and TL)
  for(i in seq_along(x)){
    list$x$sep[[i]] <- scale(x[[i]])
  }
  for(i in seq_along(y)){
    list$y$sep[[i]] <- scale(y[[i]])
  }
  # average (MTL) or concatenate (TL) data
  if(mode=="multiple"){
    list$y$com <- rowSums(sapply(X=y,FUN=scale))
    list$x$com <- x[[1]]
  } else if(mode=="transfer"){
    list$y$com <- unlist(lapply(X=y,FUN=scale))
    list$x$com <- do.call(what="rbind",args=x)
  }
  return(list)
}

init.coef <- function(list,alpha=0.95,lambda.sep=NULL,lambda.com=NULL){
  if(is.null(lambda.sep)!=is.null(lambda.com)){stop()}
  cv <- is.null(lambda.sep) & is.null(lambda.com)
  q <- length(list$x$sep)
  glm.sep <- coef.sep <- list()
  if(!is.na(alpha) & cv){
    lambda.sep <- numeric() 
  }
  # separate models
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
  # common model
  if(is.na(alpha)){
    n <- sapply(X=list$x$sep,FUN=nrow)
    #stop() # check dimensions!
    coef.com <- sapply(X=coef.sep,function(x) x) %*% n/sum(n) #
  } else if(cv){
    glm.com <- glmnet::cv.glmnet(x=list$x$com,y=list$y$com,alpha=alpha)
    coef.com <- stats::coef(glm.com,s="lambda.min")[-1]
    lambda.com <- glm.com$lambda.min
  } else {
    glm.com <- glmnet::glmnet(x=list$x$com,y=list$y$com,alpha=alpha)
    coef.com <- stats::coef(glm.com,s=lambda.com)[-1]
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

devel <- function(x,y,family="gaussian",alpha=1,nfolds=10){
  if(any(family!="gaussian")){stop("not implemented")}
  if(alpha!=1){stop("not implemented")}
  
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
    family <- rep(family,times=q)
  }
  
  list.ext <- fuse(x=x,y=y,mode=mode)
  init.ext <- init.coef(list=list.ext)

  #ncand <- 11
  #prop <- seq(from=0,to=1,length.out=ncand)
  exp <- c(0,0.5,1,2,10) # 
  exp <- seq(from=0,to=1,by=0.2)
  grid <- expand.grid(sep=exp,com=exp)
  ncand <- nrow(grid)
  model.ext <- list()
  for(i in seq_len(length.out=q)){
    model.ext[[i]] <- list()
    for(j in seq_len(ncand)){
      pf <- penfac(sep=init.ext$sep[[i]],com=init.ext$com,exp.sep=grid$sep[j],exp.com=grid$com[j])
      if(all(is.infinite(pf))){
        model.ext[[i]][[j]] <- glmnet::glmnet(x=cbind(x[[i]],-x[[i]]),y=y[[i]],family=family[i],lambda=99e99)
      } else {
        model.ext[[i]][[j]] <- glmnet::glmnet(x=cbind(x[[i]],-x[[i]]),y=y[[i]],family=family[i],lower.limits=0,penalty.factor=pf)
      }
    }
  }
  
  # create empty matrices
  pred <- list()
  for(j in seq_len(q)){
    pred[[j]] <- list()
    for(k in seq_len(ncand)){
      pred[[j]][[k]] <- matrix(NA,nrow=n[j],ncol=length(model.ext[[j]][[k]]$lambda))
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
    init.int <- init.coef(list=list.int,lambda.sep=init.ext$lambda.sep,lambda.com=init.ext$lambda.com)

    for(j in seq_len(q)){
      cond <- foldid[[j]]==i
      for(k in seq_len(ncand)){
        pf <- penfac(sep=init.int$sep[[j]],com=init.int$com,exp.sep=grid$sep[k],exp.com=grid$com[k]) # prop=prop[k],
        if(all(is.infinite(pf))){
          model.int <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]])[!cond,],y=y[[j]][!cond],family=family[j],lambda=99e99)
        } else {
          model.int <- glmnet::glmnet(x=cbind(x[[j]],-x[[j]])[!cond,],y=y[[j]][!cond],family=family[j],lower.limits=0,penalty.factor=pf)
        }
        pred[[j]][[k]][cond,] <- stats::predict(object=model.int,newx=cbind(x[[j]],-x[[j]])[cond,],s=model.ext[[j]][[k]]$lambda)
      }
    }
  }
  
  # calculate MSE
  mse <- list()
  id.grid <- lambda.min <- numeric()
  for(j in seq_len(q)){
    mse[[j]] <- list()
    for(k in seq_len(ncand)){
      mse[[j]][[k]] <- apply(X=pred[[j]][[k]],MARGIN=2,FUN=function(x) mean((x-y[[j]])^2))
    }
    id.grid[j] <- which.min(sapply(mse[[j]],min))
    lambda.min[j] <- model.ext[[j]][[id.grid[j]]]$lambda[which.min(mse[[j]][[id.grid[j]]])]
    graphics::plot(sapply(X=mse[[j]],FUN=min),type="o")
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

#object <- devel(x=x,y=y,family="gaussian")
#y_hat <- predict(object=object,newx=x)
