doParallel::registerDoParallel(ncores)
doRNG::registerDoRNG(20231118)

pfilterLikes <- data.frame(
  "ll" = rep(0, search$N_EVAL*search$N_INTERVAL),
  "which" = rep(1:search$N_INTERVAL, each = search$N_EVAL)
)

print(paste("Getting likelihoods for each value of", profile_para, "in profiling")) 
t_profile_bpf <- system.time(
  pfilterLikes$ll <- foreach(  # Get log-likelihood for each unit and set of parameters, N_EVAL times each
    i = 1:(search$N_EVAL*search$N_INTERVAL), .combine = rbind, .packages = 'spatPomp'
  ) %dopar% {
    p3 <- guesses_all[(i-1) %/% search$N_EVAL + 1, ]

    coef(model_spat) <- p3
    logLik(
      bpfilter(
        model_spat, Np = search$NP_EVAL,
        block_size = 1
      )
    )
  }
)
print(paste("Finished getting likelihoods for each value of", profile_para, "in profiling with time:"))
print(t_profile_bpf)
print(date())

# Group by model parameter set.
profile_logLik <- pfilterLikes |>
  dplyr::group_by(which) |>
  dplyr::summarize(logLik = logmeanexp(ll),
            se = logmeanexp(ll, se = TRUE)[2])

search_results <- list()
search_results$logLiks <- profile_logLik
search_results$params <- guesses_all
search_results$bpf_time <- t_profile_bpf

print(paste("Log-likelihoods for each value of", profile_para, "in profile search:"))
print(profile_logLik)

rm(profile_logLik, guesses, pfilterLikes, t_profile_bpf, guesses_all)
gc()
doParallel::stopImplicitCluster()