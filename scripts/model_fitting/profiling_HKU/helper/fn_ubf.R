ubfLikes <- data.frame(
  "ll" = rep(0, search$N_INTERVAL)
)
print(paste("Getting likelihoods for each value of", profile_para, "in profiling")) 
t_profile_ubf <- system.time(
  ubfLikes$ll <- vapply(1:search$N_INTERVAL, function(i){
    print(paste("Getting likelihood of ", profile_para, " for interval", i))
    print(date())
    doParallel::registerDoParallel(ncores)
    doRNG::registerDoRNG(20231118)

    coef(model_spat) <- guesses_all[i, ]
    log_lik <- logLik(abf(model_spat, Np=1, Nrep=search$NREPS, nbhd=fn_nbhd)) # ubf is a special case when abf with Np=1

    doParallel::stopImplicitCluster()
    return(log_lik)
  }, numeric(1))
)
print(paste("Finished getting likelihoods for each value of", profile_para, "in profiling with time:"))
print(t_profile_ubf)
print(date())

search_results <- list()
search_results$logLiks <- ubfLikes
search_results$params <- guesses_all
search_results$ubf_time <- t_profile_ubf

print(paste("Log-likelihoods for each value of", profile_para, "in profile search:"))
print(paste0(1:nrow(ubfLikes), ": ", ubfLikes$ll))

rm(guesses, ubfLikes, t_profile_ubf, guesses_all)
gc()
