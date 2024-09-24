print("Starting local search")
print(date())
t_iubf_local <- system.time(
  local_iubf <- iubf(
      object=model_spat,
      Nubf=search$Nubf,
      Nrep_per_param=search$Nrep_per_param, 
      Nparam=search$Nparam,
      nbhd=search$nbhd,
      prop=search$prop,
      sharedParNames = shared_param_names_est,
      spat_regression = search$SPAT_REGRESSION,
      rw.sd=rw.sd_s,
      cooling.type=search$cooling.type,
      cooling.fraction.50=search$cooling.fraction.50,
      verbose = TRUE
    )
)
doParallel::stopImplicitCluster()

print("Finished local search with time:")
print(t_iubf_local)
saveRDS(local_iubf, paste0(output_dir, "/local_iubf_", partition, "_", task_id, ".rds"))
local_iubf <- readRDS(paste0(output_dir, "/local_iubf_", partition, "_", task_id, ".rds"))
print(paste0("Size of local_iubf: ", format(object.size(local_iubf), units = "Gb"), " bytes"))
print(date())

iubf_params <- coef(local_iubf)
saveRDS(iubf_params, paste0(output_dir, "/local_iubf_params_", partition, "_", task_id, ".rds"))
print(paste0("Size of iubf_params: ", format(object.size(iubf_params), units = "Mb"), " bytes"))
rm(rw.sd_s, local_iubf)
gc()

print("Getting likelihoods in local search")
print(date())
doParallel::registerDoParallel(ncores)
doRNG::registerDoRNG(20231118)
coef(model_spat) <- iubf_params
t_local_ubf <- system.time(
  ubfLikes <-
    logLik(abf(model_spat, Np=1, Nrep=search$Nrep_per_param, nbhd=search$nbhd)) # ubf is a special case when abf with Np=1
)
print("Finished getting likelihoods in local search with time:")
print(t_local_ubf)
print(date())
doParallel::stopImplicitCluster()

search_results <- list()
search_results$logLiks <- ubfLikes
search_results$params <- iubf_params
search_results$iubf_time <- t_iubf_local

print("Log-likelihoods for each parameter set in local search:")
print(ubfLikes)