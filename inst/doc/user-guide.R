## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(offlineChange)
#  # Data
#  N <- 1000
#  N1 <- floor(0.1*N)
#  N2 <- floor(0.3*N)
#  a1 <- c(0.8, -0.3); c1 <- 0
#  a2 <- c(-0.5, 0.1); c2 <- 0
#  a3 <- c(0.5, -0.5); c3 <- 0
#  y <- rep(0,N)
#  L<-2
#  y[1:L] <- rnorm(L)
#  for (n in (L+1):N){
#    if (n <= N1) {
#      y[n] <- y[(n-1):(n-L)] %*% a1 + c1 + rnorm(1)
#    } else if (n <= (N1+N2)) {
#      y[n] <- y[(n-1):(n-L)] %*% a2 + c2 + rnorm(1)
#    }
#    else {
#      y[n] <- y[(n-1):(n-L)] %*% a3 + c3 + rnorm(1)
#    }
#  }
#  result <- MultiWindow(y,
#                        window_list = c(100,50,20,10,5),
#                        point_max   = 5)

## ----eval=FALSE----------------------------------------------------------
#  result$n_peak_range
#  result$peak_range

## ----eval=FALSE----------------------------------------------------------
#  result <- MultiWindow(y,
#                        window_list = c(100, 50, 20, 10, 5),
#                        point_max   = 5,
#                        prior_range = NULL,
#                        get_mle     = GetMle,
#                        penalty     = c("bic","aic"),
#                        seg_min     = 1,
#                        num_init    = NULL,
#                        tolerance   = 1,
#                        cpp         = TRUE,
#                        ret_score   = FALSE
#                        )

## ----eval=FALSE----------------------------------------------------------
#  result <- MultiWindow(y,
#                        window_list = c(100,50,20,10,5),
#                        point_max   = 5)
#  
#  RangeToPoint(y,
#               n_peak_range = result$n_peak_range,
#               peak_range   = result$peak_range)

## ----eval=FALSE----------------------------------------------------------
#  result <- MultiWindow(y,
#                        window_list = c(100,50,20,10,5),
#                        point_max   = 5)
#  
#  ChangePointsPlot(y,result,main="plot of change points")

## ----eval=FALSE----------------------------------------------------------
#  result <- MultiWindow(y,
#                        window_list = c(100,50,20,10,5),
#                        point_max   = 5,
#                        ret_score=TRUE)
#  
#  ScorePlot(result, main="score plot")

## ----eval=FALSE----------------------------------------------------------
#  install.packages('MHadaptive')
#  library('MHadaptive')
#  gen_EFdata <- function(N){
#    q <- 4
#    Ns <- c(floor(0.2*N), floor(0.3*N), N-floor(0.2*N)-floor(0.3*N))
#    thetas <- c(0.2,0.6,1)
#    X <- c()
#    for (k in 1:3){
#      loglik <- function( x ){
#        res <- - thetas[k] * abs(x)^q
#      }
#      mh <- Metro_Hastings(loglik,
#                           rgamma(n = 1, shape = 3, rate = 1),
#                           prop_sigma = NULL,
#                           par_names = NULL,
#                           iterations = 11000,
#                           burn_in = 1000,
#                           adapt_par = c(100, 20, 0.5, 0.75),
#                           quiet = TRUE)
#      mh0 <-  mcmc_thin(mh, thin = 20)
#      X <- c(X, mh0$trace[1:Ns[k]])
#    }
#    return(X)
#  }
#  GetHle=function(x,window_size) {
#    q <- 4
#    N <- length(x)
#    n_window <- ceiling(N/window_size)
#    x_transformed <- rep(0,n_window)
#    for (n in 1:n_window) {
#      #get estimated coefficients
#      xx <- x[(1+(n-1)*window_size):min(n*window_size,N)]
#      x_transformed[n] <- q * (q-1) * sum(abs(xx)^(q-2)) / sum(q^2 * abs(xx)^(2*q-2))
#    }
#    return(x_transformed)
#  }
#  
#  N <- 1000
#  X <- gen_EFdata(N)
#  result <- MultiWindow(X,
#                        window_list=c(100,80,50,30, 20),
#                        get_mle=GetHle,
#                        point_max=3,
#                        seg_min=1,
#                        tolerance=1)
#  

## ----eval=FALSE----------------------------------------------------------
#  toy_penatly=function(num_each, wgss, D, K, N) {
#    penalty <- sum(num_each * log(wgss/num_each)) + 2 * D * K * log(N)
#    return(penalty)
#  }
#  result <- MultiWindow(y,
#                        window_list=c(100,80,50,30, 20),
#                        penalty="toy_penalty",
#                        point_max=3,
#                        seg_min=1,
#                        tolerance=1)

## ----eval=FALSE----------------------------------------------------------
#  result <- MultiWindow(y,
#                        window_list = c(100,50,20,10,5),
#                        prior_range = list(c(30,200),c(220,400)))

## ----eval=FALSE----------------------------------------------------------
#  result <- MultiWindow(y,
#                        window_list = c(100, 50, 20, 10, 5),
#                        prior_range = list(c(30,200), c(220,400)),
#                        get_mle     = GetMle,
#                        tolerance   = 1)

## ----eval=FALSE----------------------------------------------------------
#  result <- MultiWindow(y,
#                        window_list = c(100,50,20,10,5),
#                        prior_range = list(c(30,200), c(220,400)))
#  
#  RangeToPoint(y,
#               n_peak_range = result$n_peak_range,
#               peak_range   = result$peak_range)

## ----eval=FALSE----------------------------------------------------------
#  # Data
#  a <- matrix(rnorm(40,mean=-1,sd=1), nrow=20, ncol=2)
#  b <- matrix(rnorm(120,mean=0,sd=1), nrow=60, ncol=2)
#  c <- matrix(rnorm(40,mean=1,sd=1), nrow=20, ncol=2)
#  x <- rbind(a,b,c)
#  result <- ChangePoints(x, point_max = 5)

## ----eval=FALSE----------------------------------------------------------
#  result <- ChangePoints(x,
#                         point_max = 5,
#                         penalty   = c("bic","aic"),
#                         seg_min   = 1,
#                         num_init  = NULL,
#                         cpp       = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  result <- OrderKmeans(x, K = 2)

## ----eval=FALSE----------------------------------------------------------
#  result <- OrderKmeans(x, K = 2, num_init=NULL)

## ----eval=FALSE----------------------------------------------------------
#  l1 <- c(15,25)
#  l2 <- c(75,100)
#  prior_range_x <- list(l1, l2)
#  result <- PriorRangeOrderKmeans(x, prior_range_x = list(l1,l2))

## ----eval=FALSE----------------------------------------------------------
#  result <- PriorRangeOrderKmeans(x, prior_range_x, num_init=NULL)

