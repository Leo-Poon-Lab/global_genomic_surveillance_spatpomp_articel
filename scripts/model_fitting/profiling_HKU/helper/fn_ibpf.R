print("Starting local search")
print(date())
t_ibpf_local <- system.time(
  foreach(
    i = 1:search$NREPS,
    .packages = c("spatPomp"),
    .combine = c
  ) %dopar% {

    ibpf_out <- ibpf(
      model_spat,
      Nbpf = search$NBPF,
      Np = search$NP,
      sharedParNames = shared_param_names_est,
      unitParNames = unit_specific_names_est,
      spat_regression = search$SPAT_REGRESSION,
      rw.sd = rw.sd_s,
      cooling.fraction.50 = search$COOLING,
      block_size = 1
      # params = r_params  Not yet implemented in spatPomp
    )
  } -> local_ibpf
)

print("Finished local search with time:")
print(t_ibpf_local)
saveRDS(local_ibpf, paste0(output_dir, "/local_ibpf_", partition, "_", n_cores, "_", task_id, ".rds"))
# local_ibpf <- readRDS(paste0(output_dir, "/local_ibpf_", partition, "_", n_cores, "_", task_id, ".rds"))
print(paste0("Size of local_ibpf: ", format(object.size(local_ibpf), units = "Gb"), " bytes"))
print(date())

pfilterLikes <- data.frame(
  "ll" = rep(0, search$NREPS_EVAL*search$NREPS),
  "which" = rep(1:search$NREPS, each = search$NREPS_EVAL)
)

ibpf_params <- t(coef(local_ibpf))
saveRDS(ibpf_params, paste0(output_dir, "/local_ibpf_params_", partition, "_", n_cores, "_", task_id, ".rds"))
print(paste0("Size of ibpf_params: ", format(object.size(ibpf_params), units = "Mb"), " bytes"))
rm(rw.sd_s, local_ibpf)
gc()
doParallel::stopImplicitCluster()

print("Getting likelihoods in local search")
print(date())
doParallel::registerDoParallel(ncores)
doRNG::registerDoRNG(20231118)
t_local_bpf <- system.time(
  pfilterLikes$ll <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
    i = 1:(search$NREPS_EVAL*search$NREPS), .combine = rbind, .packages = 'spatPomp'
  ) %dopar% {
    p <- ibpf_params[(i-1) %/% search$NREPS_EVAL + 1, ]
    coef(model_spat) <- p
    logLik(
      bpfilter(
        model_spat, Np = search$NP_EVAL,
        block_size = 1
      )
    )
  }
)
print("Finished getting likelihoods in local search with time:")
print(t_local_bpf)
print(date())
doParallel::stopImplicitCluster()

# Group by model parameter set.
ibpf_logLik <- pfilterLikes |>
  dplyr::group_by(which) |>
  dplyr::summarize(logLik = logmeanexp(ll),
                    se = logmeanexp(ll, se = TRUE)[2])

search_results <- list()
search_results$logLiks <- ibpf_logLik
search_results$params <- ibpf_params
search_results$ibpf_time <- t_ibpf_local

print("Log-likelihoods for each parameter set in local search:")
print(ibpf_logLik)