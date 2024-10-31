

#--- logit and sigmoid transforms ---

n <- 100
eps <- 1e-06

testthat::test_that("sigmoid function returns values in unit interval",{
  eta <- stats::rnorm(n=n)
  prob <- sigmoid(eta)
  testthat::expect_true(all(prob>=0 & prob<=1))
})

testthat::test_that("sigmoid followed by logit returns identity",{
  eta <- stats::rnorm(n=n)
  prob <- sigmoid(eta)
  eta_new <- logit(prob)
  testthat::expect_true(all(abs(eta-eta_new)<=eps))
})

testthat::test_that("logit followed by sigmoid returns identity",{
  prob <- stats::runif(n=n)
  eta <- logit(prob)
  prob_new <- sigmoid(eta)
  testthat::expect_true(all(abs(prob-prob_new)<=eps))
})

#--- mean and link functions ---

testthat::test_that("Gaussian mean and link functions return identity",{
  eta <- stats::runif(n=n)
  mu <- mean.function(eta=eta,family="gaussian")
  eta_new <- link.function(mu=mu,family="gaussian")
  testthat::expect_true(all(eta==mu) & all(eta==eta_new))
})

testthat::test_that("binomial mean followed by link function returns identity",{
  eta <- stats::rnorm(n=n)
  mu <- mean.function(eta=eta,family="binomial")
  eta_new <- link.function(mu=mu,family="binomial")
  testthat::expect_true(all(abs(eta-eta_new)<=eps))
})

testthat::test_that("binomial link followed by mean function returns identity",{
  mu <- stats::runif(n=n)
  eta <- link.function(mu=mu,family="binomial")
  mu_new <- mean.function(eta=eta,family="binomial")
  testthat::expect_true(all(abs(mu-mu_new)<=eps))
})

#--- Gaussian/binomial deviance ---

testthat::test_that("Gaussian deviance equals zero",{
  y <- y_hat <- stats::rnorm(n=100)
  metric <- calc.metric(y=y,y_hat=y_hat,family="gaussian")
  testthat::expect_true(metric==0)
})

testthat::test_that("binomial deviance equals zero",{
  y <- y_hat <- stats::rbinom(n=n,size=1,prob=0.5)
  metric <- calc.metric(y=y,y_hat=y_hat,family="binomial")
  testthat::expect_true(metric==0)
})

testthat::test_that("Gaussian deviance increases",{
  signal <- stats::rnorm(n=n)
  noise <- stats::rnorm(n=n)
  weight <- seq(from=0,to=1,by=0.1)
  y <- signal
  metric <- rep(x=NA,times=length(weight))
  for(i in seq_along(weight)){
    y_hat <- (1-weight[i])*signal + weight[i]*noise
    metric[i] <- calc.metric(y=y,y_hat=y_hat,family="gaussian")
  }
  testthat::expect_true(all(diff(metric)>0))
})

testthat::test_that("binomial deviance increases",{
  signal <- stats::rnorm(n=n)
  noise <- stats::rnorm(n=n)
  weight <- seq(from=0,to=1,by=0.1)
  y <- round(sigmoid(signal))
  metric <- rep(x=NA,times=length(weight))
  for(i in seq_along(weight)){
    y_hat <- sigmoid((1-weight[i])*signal + weight[i]*noise)
    metric[i] <- calc.metric(y=y,y_hat=y_hat,family="binomial")
  }
  testthat::expect_true(all(diff(metric)>0))
})

#--- fold identifiers ---

testthat::test_that("folds for multi-task learning match target",{
  q <- 3
  family <- "binomial"
  y <- sim.data.multiple(family=family,q=q,n0=100)$y_train
  fold <- make.folds.multi(y,family=family)
  code <- apply(X=y,MARGIN=1,FUN=function(x) paste0(x,collapse=""))
  for(i in seq_len(q)){
    table <- table(y[,i],fold)
    colSums <- colSums(table)
    rowSums <- rowSums(table)
    testthat::expect_true(all(rowSums==table(y[,i])))
  }
  # Think about verifying balance.
})

testthat::test_that("folds for transfer learning are balanced",{
  q <- 5
  family <- "binomial"
  y <- sim.data.transfer(family=family,q=q,n0=c(5,10,15,20,50))$y_train
  fold <- make.folds.trans(y,family=family)
  for(i in seq_len(q)){
    table <- table(y[[i]],fold[[i]])
    colSums <- colSums(table)
    rowSums <- rowSums(table)
    testthat::expect_true(all(rowSums==table(y[[i]])))
    testthat::expect_true(max(colSums)-min(colSums)<=2)
  }
})

#--- fuse data ---

testthat::test_that("fusing target works",{
  q <- 5
  data <- sim.data.multiple(q=q)
  split <- apply(X=data$y_train,MARGIN=2,FUN=function(x) x,simplify=FALSE)
  fuse <- fuse.data(x=NULL,y=split)
  testthat::expect_true(all(data$y_train==fuse$y))
})

#--- sensitivity, specificity, precision ---

truth <- sample(x=c(-1,0,1),size=n,replace=TRUE)

testthat::test_that("metrics are above zero and below one if estim=random",{
  estim <- sample(truth)
  metric <- count_vector(truth=truth,estim=estim)
  testthat::expect_true(all(metric>0 & metric<1))
})

testthat::test_that("sens=spec=prec=1 if estim=truth",{
  estim <- truth
  metric <- count_vector(truth=truth,estim=estim)
  testthat::expect_true(all(metric==1))
  metric <- count_vector(truth=abs(truth),estim=abs(estim))
  testthat::expect_true(all(metric==1))
})

testthat::test_that("spec=0 if estim=1",{
  estim <- rep(x=1,times=n)
  metric <- count_vector(truth=truth,estim=estim)
  testthat::expect_true(!metric["sensitivity"] %in% c(0,1) & metric["specificity"]==0 & !metric["precision"] %in% c(0,1))
  metric <- count_vector(truth=abs(truth),estim=abs(estim))
  testthat::expect_true(metric["sensitivity"]==1 & metric["specificity"]==0 & !metric["precision"] %in% c(0,1))
})

testthat::test_that("sens=0 and prec=1 if estim=0",{
  estim <- rep(x=0,times=n)
  metric <- count_vector(truth=truth,estim=estim)
  testthat::expect_true(metric["sensitivity"]==0 & metric["specificity"]==1 & is.na(metric["precision"]))
  metric <- count_vector(truth=abs(truth),estim=abs(estim))
  testthat::expect_true(metric["sensitivity"]==0 & metric["specificity"]==1 & is.na(metric["precision"]))
})

testthat::test_that("sens=prec=0 and spec=1 if estim=-truth",{
  estim <- -truth
  metric <- count_vector(truth=truth,estim=estim)
  testthat::expect_true(all(metric[c("sensitivity","precision")]==0) & metric["specificity"]==1)
})

testthat::test_that("count_matrix agrees with count_vector",{
  p <- 5
  truth <- sample(x=c(-1,0,1),size=n,replace=TRUE)
  estim <- matrix(data=sample(x=c(-1,0,1),size=n*p,replace=TRUE),nrow=n,ncol=p)
  metric0 <- apply(X=estim,MARGIN=2,FUN=function(x) count_vector(truth=truth,estim=x))  
  metric1 <- count_matrix(truth=matrix(data=truth,nrow=n,ncol=p),estim=estim)
  testthat::expect_true(all(metric0==metric1))
})

#--- construct penalty factors ---

w_int <- stats::runif(n=n)
w_ext <- stats::runif(n=n)

testthat::test_that("w_tot=w_int+1 if v_int=1 and v_ext=0",{
  w_tot <- 1/construct_pf(w_int=w_int,w_ext=w_ext,v_int=1,v_ext=0,type="exp")
  testthat::expect_true(all(abs(w_tot-(w_int+1))<eps))
})

testthat::test_that("w_tot=w_ext+1 if v_int=1 and v_ext=0",{
  w_tot <- 1/construct_pf(w_int=w_int,w_ext=w_ext,v_int=0,v_ext=1,type="exp")
  testthat::expect_true(all(abs(w_tot-(w_ext+1))<eps))
})

testthat::test_that("w_tot=w_int+w_ext if v_int=1 and v_ext=1",{
  w_tot <- 1/construct_pf(w_int=w_int,w_ext=w_ext,v_int=1,v_ext=1,type="exp")
  testthat::expect_true(all(abs(w_tot-(w_int+w_ext)<eps)))
})

#--- ---


  