###################################################################
########################## SepBART_blk ############################
###################################################################
##### Instead of calculating the inverse of an N by N matrix, we generate multiple 100 by 100 matrices and solve them.
###################################################################
########################### Arguments #############################
###################################################################
##### Y: Outcome
##### X: n by p covariate matrix
##### W: q dimensional exposure matrix
##### Xtest: Testing locations (not be used in data analysis)
##### W1, W0: Two exposure levels under which we want to compare the outcomes
##### W0_VIM: Exposure levels under which MTE-VIMs are calculated
##### n_blocks: The number of blocks which we want to decompose an n by n large matrix into.
##### nMCMC: The total number of MCMC iterations
##### BurnIn_portion: The proportion of iterations to be discarded
##### stepsize: The step size after the burn-in period.
###################################################################

SepBART <- function(Y,X,W, Xtest = NULL,
                        W1_quantile = 0.75,
                        W0_quantile = 0.25,
                        W0_quantile_VIM = c(0.25,0.5,0.75),
                        n_blocks = NULL,
                        Mod_ind = NULL,
                        Group_ind = NULL,
                        nMCMC=2000,BurnIn_portion=0.2,stepsize=5){
  start.time.MCMC <- Sys.time()

  num_obs <- length(Y) ## Sample size
  p <- ncol(X) ## The number of covariates
  q <- ncol(W) ## The number of exposures
  BurnIn <- nMCMC*BurnIn_portion ## The number of MCMC iterations to be discarded.
  save_index <- seq(BurnIn+1,nMCMC,stepsize) ## A sequence of iteration indices to be stored.
  nsave <- length(save_index)
  v <- length(W0_quantile_VIM)

  ## If n_blocks is not given,
  ## find the total number of blocks necessary to partition into 100 by 100 sub-matrices.
  if (is.null(n_blocks)) n_blocks <- floor(num_obs/100)

  if (is.null(Mod_ind)) Mod_ind <- 1:p

  X_names <- colnames(X)
  W_names <- colnames(W)
  Mod_names <- colnames(X)[Mod_ind]

  ## We will normalize the covariates and exposures to fall between 0 and 1 using quantiles.
  ## Define the quantile function.
  trank <- function(x) {
    x_unique <- unique(x)
    x_ranks <- rank(x_unique, ties.method = "max")
    tx <- x_ranks[match(x,x_unique)] - 1

    tx <- tx / length(unique(tx))
    tx <- tx / max(tx)

    return(tx)
  }
  quantile_normalize_bart <- function(X) {
    apply(X = X, MARGIN = 2, trank)
  }

  ## We scale the outcome and store the original standard deviation.
  sd_Y <- sd(Y)
  Y_orig <- Y
  Y <- Y/sd_Y

  ## Save the quantiles of exposures and covariates
  X_quantiles_orig <- list()
  W_quantiles_orig <- list()
  for (j in 1:p){
    X_quantiles_orig[[j]] <- quantile(X[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
  }
  for (j in 1:q){
    W_quantiles_orig[[j]] <- quantile(W[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
  }

  ## default W1 and W0
  W1 = apply(W, 2, quantile, W1_quantile, na.rm = TRUE)
  W0 = apply(W, 2, quantile, W0_quantile, na.rm = TRUE)
  print(cat("W0:",W0))
  print(cat("W1:",W1))


  W0_VIM = apply(W, 2, quantile, W0_quantile_VIM, na.rm = TRUE)
  cat("W0_VIM:")
  print(W0_VIM)

  ## Calculate propensity and choose the upper half
  temp <- data.frame(W=W, X=X)
  fit <- lm(W ~ X + I(X^2) + I(X^3), data=temp)
  mean <- fit$fitted.values
  res_cov <- (t(fit$residuals)%*%fit$residuals) / (num_obs-p)

  ps1 <- ps0 <- ps_min <- rep(NA, num_obs)
  W1_window <- apply(W,2, quantile, c(0.65,0.85), na.rm = T)
  W0_window <- apply(W,2, quantile, c(0.15,0.35), na.rm = T)
  for (i in 1:num_obs){
    ps1[i] <- pmvnorm(lower=W1_window[1,],upper = W1_window[2,], mean=mean[i,], sigma = res_cov)
    ps0[i] <- pmvnorm(lower=W0_window[1,],upper = W0_window[2,], mean=mean[i,], sigma = res_cov)
  }
  ps_min <- pmin(ps1,ps0)
  upper_half <- ps_min > median(ps_min)

  if (is.null(Xtest)==FALSE){
    norm_X_Xtest <- quantile_normalize_bart(rbind(X,Xtest))
    X <- norm_X_Xtest[1:num_obs,]
    Xtest <- norm_X_Xtest[-(1:num_obs),]
  }else{
    norm_X <- quantile_normalize_bart(X)
    X <- norm_X[1:num_obs,]
    norm_W <- quantile_normalize_bart(W)
    W <- norm_W
    W1 = apply(norm_W, 2, quantile, W1_quantile, na.rm = TRUE)
    W0 = apply(norm_W, 2, quantile, W0_quantile, na.rm = TRUE)
    W0_VIM = apply(norm_W, 2, quantile, W0_quantile_VIM, na.rm = TRUE)
  }

  ## W quantiles after normalizing. This will be used to calculate E(Y(W)-Y(W0))
  X_sample_size <- min(500, num_obs)
  W_quantiles <- list()
  for (j in 1:q){
    W_quantiles[[j]] <- quantile(W[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
  }
  W_quantiles_length <- length(W_quantiles[[j]])
  W1_rep <- t(W1)[rep(1,num_obs),] ## n by q matrix with each row of W1. Use this to calculate h(X,W1)
  W0_rep <- t(W0)[rep(1,num_obs),] ## n by q matrix with each row of W0. Use this to calculate h(X,W0)
  if (is.null(Xtest)==FALSE){
    W1_rep_test <- t(W1)[rep(1,dim(Xtest)[1]),]
    W0_rep_test <- t(W0)[rep(1,dim(Xtest)[1]),]
  }


  ## Our model is
  ## E(Y|X) = f(X) + G(W) + h(X,W) where h(X,W) = sum{h_j(X_j,W)} where X_j denotes the j-th covariate.
  ## To make the model identifiable, we define X* and W* as the reference vectors
  ## so that f(X*) = G(W*) = 0 and h(X*, . ) = h( . , W*) = 0.
  ## X* and W* can be considered as the baseline of the model.
  ## We define X* and W* to be the sample (column) means of X and W.
  X_star <- apply(X,2,function(x) sort(x)[floor(num_obs/2)])
  W_star <- apply(W,2,function(x) sort(x)[floor(num_obs/2)])

  ## The interaction term h(.,.) equals 0 for any rows of three matrices below.
  ## n by (p+q) matrix with the i-th row of (W_i, X*) where W_i denotes the i-th observation's exposure vector.
  obs_grid_x0 <- 0
  for (i in 1:num_obs){
    obs_grid_x0 <- rbind(obs_grid_x0, c(W[i,],X_star))
  }
  obs_grid_x0 <- obs_grid_x0[-1,]
  ## Similarly define the n by (p+q) matrix with the i-th row of (W*, X_i).
  obs_grid_w0 <- 0
  for (i in 1:num_obs){
    obs_grid_w0 <- rbind(obs_grid_w0, c(W_star,X[i,]))
  }
  obs_grid_w0 <- obs_grid_w0[-1,]
  ## n by (p+q) matrix with each row of (W*,X*)
  obs_grid_00 <- matrix(rep(c(W_star,X_star),num_obs), nrow = num_obs, byrow = T)

  block_index <- sample(rep(1:n_blocks, len=num_obs))
  X_sampled <- list()
  W_sampled <- list()
  X_sampled_grid <- list()
  W_sampled_grid <- list()
  W_star_sampled <- list()
  X_star_sampled <- list()
  W0_sampled <- W0_grid <- list()
  for (block in 1:n_blocks){
    sample_index <- which(block_index==block)
    N <- sum(block_index==block)
    X_sampled[[block]] <- X[sample_index,]
    W_sampled[[block]] <- W[sample_index,]
    X_sampled_grid[[block]] <- X_sampled[[block]][rep(1:N,N),]
    W_sampled_grid[[block]] <- W_sampled[[block]][rep(1:N,each=N),]
    W_star_sampled[[block]] <- t(W_star)[rep(1, N),]
    X_star_sampled[[block]] <- t(X_star)[rep(1, N),]
    W0_grid[[block]] <- W0_sampled[[block]] <- vector("list",v)
    for (vv in 1:v){
      W0_grid[[block]][[vv]] <- t(W0_VIM[vv,])[rep(1,N^2),]
      W0_sampled[[block]][[vv]] <- t(W0_VIM[vv,])[rep(1, N),]
    }
  }


  ## Set the starting values
  mu_hat_c_prev <- rep(0,num_obs)
  mu_hat_X_prev <- rep(0,num_obs)
  mu_hat_W_prev <- rep(0,num_obs)
  mu_hat_intact_prev <- list()
  mu_hat_intact_new <- list()
  for (j in Mod_ind){
    mu_hat_intact_prev[[j]] <- rep(0,num_obs)
  }

  ## We will also calculate
  ## CATE_hat: E(Y(W1)-Y(W0)|X)
  ## psi_j = (Heterogeneity without X_j)/(Total heterogeneity): Variable importance measures
  if (is.null(Xtest)==FALSE){
    Xtest <- as.matrix(Xtest)
    CATE_test_hat <- E_Y_W0_test <- E_Y_W1_test <- matrix(NA, nrow=nsave, ncol=dim(Xtest)[1])
  }else{CATE_test_hat <- E_Y_W0_test <- E_Y_W1_test <- NA}
  CATE_hat <- E_Y_W1 <- E_Y_W0 <- matrix(NA, nrow=nsave, ncol=num_obs)
  Hetero_total <- matrix(NA, nsave, v)
  Hetero_wo <- Hetero_wo_group <- list()
  psi_hat <- psi_hat_group <- list()
  CATE_h_hat <- list()
  for (j in Mod_ind){
    CATE_h_hat[[j]] <- matrix(NA, nrow=nsave, ncol=num_obs)
    Hetero_wo[[j]] <- matrix(NA, nsave, v)
    psi_hat[[j]] <- matrix(NA, nsave, v)
  }
  if (!is.null(Group_ind)){
    for (gg in 1:length(Group_ind)){
      Hetero_wo_group[[gg]] <- matrix(NA, nsave, v)
      psi_hat_group[[gg]] <- matrix(NA, nsave, v)
    }
  }


  ## Estimate the conditional distribution X_j|X_{-j} using the normality assumption
  ## This will be used to calculate the heterogeneity without the j-th covariate.
  normal_weight <- normal_weight_group <- list()
  for (j in Mod_ind){
    normal_cd <- lm(X[,j]~X[,-j])
    normal_weight[[j]] <- list()
    for (block in 1:n_blocks){
      N <- sum(block_index==block)
      normal_weight[[j]][[block]] <- matrix(NA, nrow = N, ncol = N)
      for (w in 1:N){
        normal_weight[[j]][[block]][w,] <- dnorm(X_sampled[[block]][,j], mean = c(1,X_sampled[[block]][w,-j]) %*% normal_cd$coefficients, sd=sigma(normal_cd))
        normal_weight[[j]][[block]][w,] <- normal_weight[[j]][[block]][w,]/sum(normal_weight[[j]][[block]][w,])
      }
    }
  }
  if(!is.null(Group_ind)){
    for (gg in 1:length(Group_ind)){
      group <- Group_ind[[gg]]
      X_group <- X[,group]
      X_rem <- X[,-group]

      normal_cd_group <- lm(X_group~X_rem)
      res_cov <- (t(normal_cd_group$residuals)%*%normal_cd_group$residuals) / (dim(X_rem)[1]-dim(X_rem)[2])
      normal_weight_group[[gg]] <- list()
      for (block in 1:n_blocks){
        N <- sum(block_index==block)
        normal_weight_group[[gg]][[block]] <- matrix(NA, nrow = N, ncol = N)
        for (w in 1:N){
          normal_weight_group[[gg]][[block]][w,] <- dmvnorm(X_sampled[[block]][,group], mean = c(1,X_sampled[[block]][w,-group]) %*% normal_cd_group$coefficients, sigma=res_cov)
          normal_weight_group[[gg]][[block]][w,] <- normal_weight_group[[gg]][[block]][w,]/sum(normal_weight_group[[gg]][[block]][w,])
        }
      }
    }
  }
  ## Start MCMC
  ## We use the Bayesian backfitting method.
  for (k in 2:nMCMC){
    ## Update the constant term "c"
    R_int_a <- 0 ## The sum of the interactions not updated for the current scan
    R_int_b <- 0 ## The sum of the interactions updated for the current scan
    for (j in Mod_ind) {R_int_a <- R_int_a + mu_hat_intact_prev[[j]]} ## All interactions are not updated.
    R <- Y - mu_hat_X_prev - mu_hat_W_prev - R_int_a - R_int_b ## Calculate the residuals w.o. the constant term.
    mu_hat_c_new <- mu_hat_c_prev + mean(R) ## Update "c" by taking the mean of residuals.

    ## Update g(W)
    R <- Y - mu_hat_c_new - mu_hat_X_prev - R_int_a - R_int_b ## Residuals w.o. g(W).
    if (k==2) { ## Make a tree when we start MCMC.
      opts_soft_W <- SoftBart::Opts(num_burn = 0, num_save = 1, num_thin = 1, update_sigma = TRUE,
                                    update_sigma_mu = FALSE) ## Options for the tree
      hypers_W <- SoftBart::Hypers(W,Y)
      forest_W <- SoftBart::MakeForest(hypers_W,opts_soft_W)
    }
    mu_hat_W_new <- t(forest_W$do_gibbs(W,R,W,1)) ## Update g(W) by do_gibbs(explanatory, response, evaluated at, # of scans)
    mu_hat_W_w0 <- forest_W$do_predict(matrix(W_star,nrow=1)) ## Calculate g(W*)

    shared_sigma <- forest_W$get_sigma() ## f, g, h will share the same residual variance.

    ## Update f(X)
    R <- Y - mu_hat_c_new - mu_hat_W_new - R_int_a - R_int_b
    if (k==2) {
      opts_soft_X <- SoftBart::Opts(num_burn = 0, num_save = 1, num_thin = 1, update_sigma = TRUE,
                                    update_sigma_mu = FALSE) ## Options for the tree
      hypers_X <- SoftBart::Hypers(X,Y, num_tree = 50, sigma_hat = shared_sigma) ## f will use the residual variance from g.
      forest_X <- SoftBart::MakeForest(hypers_X,opts_soft_X)
    }else{
      forest_X$set_sigma(shared_sigma) ## f will use the residual variance from g.
    }
    mu_hat_X_new <- t(forest_X$do_gibbs(X,R,X,1)) ## Update f(X)
    mu_hat_X_x0 <- forest_X$do_predict(matrix(X_star,nrow=1)) ## Calculate f(X*)

    ## Update h(X,W) = sum{h_j(X_j,W)}
    if (k==2) {
      opts_mod <- ModBart::Opts(num_burn = 0, num_save = 1, num_thin = 1, update_sigma = FALSE)
      hypers <- ModBart::Hypers(W,Y, sigma_hat = shared_sigma)

      hypers$length_scale <- 4 / pi
      mean_ell_sq <- hypers$length_scale^2
      hypers$shape_length_scale <- 1
      hypers$rate_length_scale <- 1 / (mean_ell_sq)

      forest <- list()
      for (j in Mod_ind) {forest[[j]] <- ModBart::MakeForest(hypers = hypers, opts = opts_mod)}
    }else {
      for (j in Mod_ind) {forest[[j]]$set_sigma(shared_sigma)}
    }
    rd_order <- sample(Mod_ind,length(Mod_ind)) ## Random order for the update

    mu_hat_intact_grid_x0 <- list()
    mu_hat_intact_grid_w0 <- list()
    mu_hat_intact_grid_00 <- list()

    temp_count = c()
    temp_sigma_mu = c()
    for (j in rd_order){
      R_int_a <- R_int_a - mu_hat_intact_prev[[j]]
      R <- Y - mu_hat_c_new - mu_hat_X_new - mu_hat_W_new - R_int_a - R_int_b
      mu_hat_intact_new[[j]] <- t(forest[[j]]$do_gibbs(W, R, X[,j], W, X[,j], 1))
      R_int_b <- R_int_b + mu_hat_intact_new[[j]]

      mu_hat_intact_grid_x0[[j]] <- forest[[j]]$predict(obs_grid_x0[,1:q],obs_grid_x0[,q+j]) ## Calculate h_j(X*_j, W)
      mu_hat_intact_grid_w0[[j]] <- forest[[j]]$predict(obs_grid_w0[,1:q],obs_grid_w0[,q+j]) ## Calculate h_j(X_j, W*)
      mu_hat_intact_grid_00[[j]] <- forest[[j]]$predict(obs_grid_00[,1:q],obs_grid_00[,q+j]) ## Calculate h_j(X*_j, W*)
      temp_count = cbind(temp_count, forest[[j]]$get_counts())
      temp_sigma_mu = cbind(temp_sigma_mu,forest[[j]]$get_params()$sigma_mu)
    }


    ## For the iterations that will not be discarded, calculate CATE E(Y(W1)-Y(W0)|X), the variable importance measures (VIM), psi_j, for each covariate.
    ## psi_j = (Heterogeneity without X_j)/(Total heterogeneity)
    ## where Total heterogeneity = E_W[Var_X{E(Y(W)-Y(W0)|X)}], the average (over W) variability (over X) of CATE with fixed W0.
    ## We similarly define "Heterogeneity without X_j" but integrating out X_j.
    ## We calculate the VIM for each sub-matrix and take the average of VIMs.
    if (k %in% save_index){
      k_index <- which(save_index==k)
      ## Calculate CATE E(Y(W1)-Y(W0)|X) and CATE_hj E(Y(W1)-Y(W0)|X_j)
      ## In our model E(Y(W1)-Y(W0)|X) = g(W1)-g(W0) + h(X,W1)-h(X,W0)
      ## Also calculate E(Y(W)-Y(W0)|X=EX) = g(W)-g(W0)

      h_W0 <-0
      h_W1 <-0
      for (j in Mod_ind) {
        temp1 <- forest[[j]]$predict(W1_rep, X[,j])
        temp0 <- forest[[j]]$predict(W0_rep, X[,j])
        h_W1 <- h_W1 + temp1
        h_W0 <- h_W0 + temp0
        CATE_h_hat[[j]][k_index,] <- (temp1 - temp0) -
          rep(forest[[j]]$predict(t(W1), X_star[j]) - forest[[j]]$predict(t(W0), X_star[j]), num_obs)
      }
      E_Y_W1[k_index,] <- mu_hat_c_new + forest_X$do_predict(X) +
        rep(forest_W$do_predict(t(W1)),num_obs) + h_W1
      E_Y_W0[k_index,] <- mu_hat_c_new + forest_X$do_predict(X) +
        rep(forest_W$do_predict(t(W0)),num_obs) + h_W0

      CATE_hat[k_index,] <- E_Y_W1[k_index,] - E_Y_W0[k_index,]
      cat("ATE = ", mean(CATE_hat[k_index,])* sd_Y, "E_Y_W0=", mean(E_Y_W0[k_index,])* sd_Y, "E_Y_W1=", mean(E_Y_W1[k_index,])* sd_Y, "\n")

      if (is.null(Xtest)==FALSE){
        CATE_test_hat_int <- 0
        for (j in 1:p) {
          CATE_test_hat_int <- CATE_test_hat_int + forest[[j]]$predict(W1_rep_test, Xtest[,j]) - forest[[j]]$predict(W0_rep_test, Xtest[,j])
        }
        CATE_test_hat[k_index,] <- c(forest_W$do_predict(t(W1))-forest_W$do_predict(t(W0))) + CATE_test_hat_int
      }

      Hetero_total_blk <- array(NA, dim=c(n_blocks,v))
      Hetero_wo_blk <- array(NA, dim=c(p, n_blocks, v))
      Hetero_wo_blk_group <- array(NA, dim=c(length(Group_ind), n_blocks, v))
      for (block in 1:n_blocks){
        ## First, calculate the total heterogeneity.
        ## E(Y(W)-Y(W0)|X) = g(W)-g(W0) + h(X,W)-h(X,W0)
        ##                 = g_tilde + h_tilde.
        ## Since g and h (also f) are not yet identifiable, we make some adjustments during the calculation.
        N <- sum(block_index==block)
        g_tilde_int1 <- 0
        g_tilde_int2 <- rep(0, v)
        for (j in Mod_ind) {
          g_tilde_int1 <- g_tilde_int1 + forest[[j]]$predict(W_sampled[[block]], X_star_sampled[[block]][,j])
          for (vv in 1:v) {
            g_tilde_int2[vv] <- g_tilde_int2[vv] + forest[[j]]$predict(t(W0_VIM[vv,]), X_star[j])

          }
        }
        g_tilde <- matrix(NA, N , v)
        for (vv in 1:v){
          g_tilde[,vv] <-
            (forest_W$do_predict(W_sampled[[block]]) + g_tilde_int1 ) -
            (rep(forest_W$do_predict(t(W0_VIM[vv,])) + g_tilde_int2[vv], N) )
        }

        h_tilde <- list()
        for (j in Mod_ind){
          temp1 <- (forest[[j]]$predict(W_sampled_grid[[block]], X_sampled_grid[[block]][,j]) -
                      rep(forest[[j]]$predict(W_sampled[[block]], X_star_sampled[[block]][,j]), each=N ) )
          h_tilde[[j]] <- matrix(NA, N^2, v)
          for (vv in 1:v){
            h_tilde[[j]][,vv] <- temp1 -
              (forest[[j]]$predict(W0_grid[[block]][[vv]], X_sampled_grid[[block]][,j]) -
                 rep(forest[[j]]$predict(W0_sampled[[block]][[vv]], X_star_sampled[[block]][,j]), each=N) )
          }

        }

        ## Take the mean of variance of CATE=(g_tilde + h_tilde).
        Hetero_vec_int <- rep(0,v)
        for (j in Mod_ind){
          Hetero_vec_int <- Hetero_vec_int + h_tilde[[j]]
        }
        Hetero_total_vec <-
          rep(g_tilde, each=N) + Hetero_vec_int # c(tau(X1,W1),tau(X2,W1),...,tau(X1,W2),...,tau(XN,WN))
        Hetero_total_mat <- array(NA, dim=c(N,N,v))
        for (vv in 1:v){
          Hetero_total_mat[,,vv] <- matrix(Hetero_total_vec[,vv], nrow = N)
        }
        Hetero_total_blk[block,] <- apply(apply(Hetero_total_mat, c(2,3), var),2,mean)

        ## Second, calculate the heterogeneity in E(Y(W)-Y(W0)|X) without X_j.
        ## We are going to calculate E{E(Y(W)-Y(W0)|X)|X_{-j}}
        ##    where the outer conditional expectation is w.r.t. the conditional distribution of X_j|X_{-j}
        ##    where X_{-j} denotes the (p-1) covariates without X_j.
        for (j in Mod_ind){
          ## Assume X_j|X_{-j}~N(.,.) where the mean and the variance are estimated from the linear regression model.
          ## Integrating out X_j only affects h_tilde[[j]].
          ## Hence, we decompose E{E(Y(W)-Y(W0)|X)|X_{-j}} into two parts.

          ## (I) the part that does not involve X_j. -> remains as it is.
          Hetero_wo_j_vec_int = Hetero_vec_int - h_tilde[[j]]
          Hetero_wo_j_vec <-
            rep(g_tilde, each=N) + Hetero_wo_j_vec_int
          Hetero_wo_j_mat_base <- array(NA, dim=c(N,N,v))
          for (vv in 1:v){
            Hetero_wo_j_mat_base[,,vv] <- matrix(Hetero_wo_j_vec[,vv], nrow = N)
          }

          ## (II) the part that does involve X_j. -> integrated out.
          h_j_tilde_mat <- array(NA, dim=c(N,N,v))
          E_h_j_tilde_mat <- array(NA, dim=c(N,N,v))
          for (vv in 1:v){
            h_j_tilde_mat[,,vv] <- matrix(h_tilde[[j]][,vv], nrow = N)
            E_h_j_tilde_mat[,,vv] <- normal_weight[[j]][[block]] %*% h_j_tilde_mat[,,vv] # Integrating out the CATE using the conditional density.
          }

          Hetero_wo_j_mat_ind <- Hetero_wo_j_mat_base + E_h_j_tilde_mat # n by n matrix
          Hetero_wo_blk[j,block,] <- apply(apply(Hetero_wo_j_mat_ind, c(2,3), var),2,mean) # Take the variance for each W value and then take the mean of n variances.
        }
        if (!is.null(Group_ind)){
          for (gg in 1:length(Group_ind)){
            ## Assume X_j|X_{-j}~N(.,.) where the mean and the variance are estimated from the linear regression model.
            ## Integrating out X_j only affects h_tilde[[j]].
            ## Hence, we decompose E{E(Y(W)-Y(W0)|X)|X_{-j}} into two parts.
            group <- Group_ind[[gg]]
            rem <- setdiff(1:p, group)
            ## (I) the part that does not involve X_j. -> remains as it is.
            h_tilde_group <- 0
            for (j in group){
              h_tilde_group <- h_tilde_group + h_tilde[[j]]
            }
            Hetero_wo_j_vec_int_group = Hetero_vec_int - h_tilde_group
            Hetero_wo_j_vec_group <-
              rep(g_tilde, each=N) + Hetero_wo_j_vec_int_group
            Hetero_wo_j_mat_base_group <- array(NA, dim=c(N,N,v))
            for (vv in 1:v){
              Hetero_wo_j_mat_base_group[,,vv] <- matrix(Hetero_wo_j_vec_group[,vv], nrow = N)
            }

            ## (II) the part that does involve X_j. -> integrated out.
            h_j_tilde_mat_group <- array(NA, dim=c(N,N,v))
            E_h_j_tilde_mat_group <- array(NA, dim=c(N,N,v))
            for (vv in 1:v){
              h_j_tilde_mat_group[,,vv] <- matrix(h_tilde_group[,vv], nrow = N)
              E_h_j_tilde_mat_group[,,vv] <- normal_weight_group[[gg]][[block]] %*% h_j_tilde_mat_group[,,vv] # Integrating out the CATE using the conditional density.
            }

            Hetero_wo_j_mat_ind_group <- Hetero_wo_j_mat_base_group + E_h_j_tilde_mat_group # n by n matrix
            Hetero_wo_blk_group[gg,block,] <- apply(apply(Hetero_wo_j_mat_ind_group, c(2,3), var),2,mean) # Take the variance for each W value and then take the mean of n variances.
          }
        }
      }
      Hetero_total[k_index,] <- apply(Hetero_total_blk,2,mean) # Take the mean of the total heterogeneity.
      if (!is.null(Group_ind)){
        for (gg in 1:length(Group_ind)){
          Hetero_wo_group[[gg]][k_index,] <- apply(Hetero_wo_blk_group[gg,,],2,mean) # Take the mean of the heterogeneity w.o X_j.
          psi_hat_group[[gg]][k_index,] <- apply(1-Hetero_wo_blk_group[gg,,]/Hetero_total_blk, 2, mean) # Take the ratio.
        }
      }
      for (j in Mod_ind){
        Hetero_wo[[j]][k_index,] <- apply(Hetero_wo_blk[j,,], 2, mean) # Take the mean of the heterogeneity w.o X_j.
        psi_hat[[j]][k_index,] <- apply(1-Hetero_wo_blk[j,,]/Hetero_total_blk, 2, mean) # Take the ratio.
      }
    }


    #### shifting ####
    ## This step is needed to make all functions identifiable
    ##    in the sense that f(X*) = G(W*) = 0 and h(X*, . ) = h( . , W*) = 0.
    sum_mu_hat_intact_grid_00 <- 0
    sum_mu_hat_intact_grid_x0 <- 0
    sum_mu_hat_intact_grid_w0 <- 0

    for (j in Mod_ind){
      sum_mu_hat_intact_grid_00 <- sum_mu_hat_intact_grid_00 + mu_hat_intact_grid_00[[j]]
      sum_mu_hat_intact_grid_x0 <- sum_mu_hat_intact_grid_x0 + mu_hat_intact_grid_x0[[j]]
      sum_mu_hat_intact_grid_w0 <- sum_mu_hat_intact_grid_w0 + mu_hat_intact_grid_w0[[j]]
    }

    mu_hat_c_new <- mu_hat_c_new + c(mu_hat_X_x0) + c(mu_hat_W_w0) + sum_mu_hat_intact_grid_00
    mu_hat_X_new <- mu_hat_X_new - c(mu_hat_X_x0) + sum_mu_hat_intact_grid_w0 - sum_mu_hat_intact_grid_00
    mu_hat_W_new <- mu_hat_W_new - c(mu_hat_W_w0) + sum_mu_hat_intact_grid_x0 - sum_mu_hat_intact_grid_00
    for (j in Mod_ind){
      mu_hat_intact_new[[j]] <- mu_hat_intact_new[[j]] - mu_hat_intact_grid_x0[[j]] - mu_hat_intact_grid_w0[[j]] + mu_hat_intact_grid_00[[j]]
    }

    mu_hat_c_prev <- mu_hat_c_new
    mu_hat_X_prev <- mu_hat_X_new
    mu_hat_W_prev <- mu_hat_W_new
    mu_hat_intact_prev[[j]] <- mu_hat_intact_new[[j]]

    if(k%%100==0) cat(k, "/", nMCMC, "Done", "\r")

  }

  print("MCMC scan done")

  ## Burning and stepping
  X_quantiles <- list()
  for (j in Mod_ind){
    X_quantiles[[j]] <- quantile(X[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
    Index_quantile <- match(X_quantiles[[j]], X[,j])
    CATE_h_hat[[j]] <- (CATE_h_hat[[j]][, Index_quantile])* sd_Y
  }

  CATE_hat <- CATE_hat * sd_Y
  E_Y_W0 <- E_Y_W0 * sd_Y
  E_Y_W1 <- E_Y_W1 * sd_Y
  if (is.null(Xtest)==FALSE){CATE_test_hat <- CATE_test_hat * sd_Y }

  ## Calculate the posterior means
  ATE_hat <- ATE_upper_hat <- rep(NA, nsave)
  for (kk in 1:nsave){
    BB_alpha <- rep(1, num_obs)
    BB_weights <- rdirichlet(1,BB_alpha)
    ATE_hat[kk] <- BB_weights %*% CATE_hat[kk,]
    BB_alpha_upper <- rep(1, sum(upper_half))
    BB_weights_upper <- rdirichlet(1,BB_alpha_upper)
    ATE_upper_hat[kk] <- BB_weights_upper %*% CATE_hat[kk,upper_half]
  }

  end.time.MCMC <- Sys.time()
  (running.time.MCMC <- end.time.MCMC-start.time.MCMC)
  print(running.time.MCMC)

  return(list(
    CATE_h_hat = CATE_h_hat,
    ATE_hat = ATE_hat,
    ATE_upper_hat = ATE_upper_hat,
    E_Y_W0 = E_Y_W0,
    E_Y_W1 = E_Y_W1,
   CATE_hat = CATE_hat,
   CATE_test_hat = CATE_test_hat,
    Hetero_total = Hetero_total,
    psi_hat = psi_hat,
    psi_hat_group = psi_hat_group,
    sd_Y=sd_Y,
    res_sd = shared_sigma,
    X_quantiles_orig = X_quantiles_orig,
    W_quantiles_orig = W_quantiles_orig,
    X_names = X_names,
    W_names = W_names,
    Mod_names = Mod_names,
    W1 = W1,
    W0 = W0,
    W0_VIM = W0_VIM,
    W1_window = W1_window,
    W0_window = W0_window,
    ps1 = ps1,
    ps0 = ps0,
    upper_half = upper_half,
    Group_ind = Group_ind
  ))

}
