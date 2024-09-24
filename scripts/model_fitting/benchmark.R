library(pomp)
set.seed(123)
model_name <- "Omicron20"
model_omicron20 <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
x <- obs(model_omicron20)
tol = 1e-300

benchdir <- paste0("results/benchmark/")
if (!dir.exists(benchdir)) dir.create(benchdir)

stew(file=paste0(benchdir,"iid.rda"),{
  iid_loglik <- function(theta, x) {
    x_int <- as.integer(x)
    
    means <- exp(theta)    
    sds <- pmax(sqrt(exp(theta)), 1)  
    
    upper_p <- pnorm(x_int + 0.5, mean = means, sd = sds)
    lower_p <- pnorm(x_int - 0.5, mean = means, sd = sds)
    
    prob_diff <- upper_p - lower_p + tol
    
    log_prob <- ifelse(x_int > tol, log(prob_diff), 0)
    
    n <- nrow(x) 
    
    weights <- rep(1, n) 
    
    mid_start <- 2 * 29 + 1
    mid_end <- n - 4 * 29
    weights[mid_start:mid_end] <- 1 / 29
    
    weights[(n - 4 * 29 + 1):(n - 3 * 29)] <- 1  
    weights[(n - 3 * 29 + 1):(n - 2 * 29)] <- 1.5  
    weights[(n - 2 * 29 + 1):(n - 1 * 29)] <- 2  
    weights[(n - 29 + 1):n] <- 1  
    
    expanded_weights <- rep(weights, ncol(x))
    
    weighted_log_prob <- log_prob * expanded_weights
    
    return(-sum(weighted_log_prob))
  }
  
  theta <- log(rowMeans(x) + tol)
  
  iid_mle <- optim(theta, iid_loglik, x = x)
  
})
-iid_mle$value


stew(file=paste0(benchdir,"ar.rda"), {
  ar_loglik <- function(theta, x) {
    x_int <- as.integer(x)
    
    means <- exp(theta) * 0.5 + 0.5 * as.numeric(cbind(exp(theta), x[, 1:(ncol(x)-1)]))  
    sds <- pmax(sqrt(means), 1)   
    
    upper_p <- pnorm(x_int + 0.5, mean = means, sd = sds)
    lower_p <- pnorm(x_int - 0.5, mean = means, sd = sds)
    
    prob_diff <- upper_p - lower_p + tol
    
    log_prob <- ifelse(x_int > tol, log(prob_diff), 0)
    
    n <- nrow(x)
    
    weights <- rep(1, n)
    
    mid_start <- 2 * 29 + 1
    mid_end <- n - 4 * 29
    weights[mid_start:mid_end] <- 1 / 29
    
    weights[(n - 4 * 29 + 1):(n - 3 * 29)] <- 1   
    weights[(n - 3 * 29 + 1):(n - 2 * 29)] <- 1.5  
    weights[(n - 2 * 29 + 1):(n - 1 * 29)] <- 2   
    weights[(n - 29 + 1):n] <- 1   
    
    expanded_weights <- rep(weights, ncol(x))
    
    weighted_log_prob <- log_prob * expanded_weights
    
    return(-sum(weighted_log_prob))
  }
  
  theta <- log(rowMeans(x) + tol)
  
  ar_mle <- optim(theta, ar_loglik, x = x)
  
})
-ar_mle$value

