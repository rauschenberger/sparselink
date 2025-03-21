
fuse <- function(x,y,mode){
  slot <- list(sep=list(),com=numeric())
  list <- list(y=slot,x=slot)
  for(i in seq_along(x)){
    list$x$sep[[i]] <- scale(x[[i]])
  }
  for(i in seq_along(y)){
    list$y$sep[[i]] <- scale(y[[i]])
  }
  if(mode=="multiple"){
    list$y$com <- rowSums(sapply(X=y,FUN=scale))
    list$x$com <- x[[1]]
  } else if(mode=="transfer"){
    list$y$com <- unlist(lapply(X=y,FUN=scale))
    list$x$com <- do.call(what="rbind",args=x)
  }
  return(list)
}


devel <- function(x,y,family,alpha.init=0.95,alpha=1,type="exp",nfolds=10){
  
  alpha.one <- alpha.init
  alpha.two <- alpha

  cat("type=",type,"\n")
  
  if(is.matrix(y) & is.matrix(x)){
    message("mode: multi-target learning")
    p <- ncol(x)
    q <- ncol(y)
    n <- rep(x=nrow(x),times=q)
    foldid <- make.folds.multi(y=y,family=family,nfolds=nfolds)
    y <- apply(y,2,function(x) x,simplify=FALSE)
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
  
  list <- fuse(x=x,y=y,mode=mode)
  glm.sep.ext <- coef.sep.ext <- list()
  for(i in seq_len(q)){
    glm.sep.ext[[i]] <- glmnet::cv.glmnet(x=list$x$sep[[i]],y=list$y$sep[[i]])
    coef.sep.ext[[i]] <- coef(glm.sep.ext[[i]],s="lambda.min")[-1]
  }
  glm.com.ext <- glmnet::cv.glmnet(x=list$x$com,y=list$y$com)
  coef.com.ext <- stats::coef(glm.com.ext,s="lambda.min")[-1]
  
  # CONTINUE HERE
  
  class(list) <- "devel"
  return(list)
}
