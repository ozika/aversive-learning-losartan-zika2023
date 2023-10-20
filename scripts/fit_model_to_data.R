### Script to estimate bayesian models and output likelihood (loo/waic) and parameter estimates

extract_visit_sub_outc <- function(df) {
  # Apply the function and create new columns
  df <- df %>%
    mutate(numbers = lapply(variable, extract_numbers)) %>%
    mutate(visit = sapply(numbers, `[`, 1),
           sub = sapply(numbers, `[`, 2),
           outcm = sapply(numbers, `[`, 3)) %>%
    select(-numbers)
  return(df)
  
}
# helper func to extract numbers
extract_numbers <- function(input_string) {
  #numbers <- as.numeric(gsub("[^0-9]", " ", input_string))
  matches <- regmatches(input_string, gregexpr("\\d+", input_string))
  # Convert extracted numbers to numeric
  numbers <- as.numeric(unlist(matches))
  return(numbers)
}

strip_ll_of_nans <- function(ll_raw) {
  nvalid<- length(ll_raw[1,!is.na(ll_raw[1,])])
  ll <- array(dim=c(dim(ll_raw)[1], nvalid))
  for (i in 1:dim(ll_raw)[1]) {
    ll[i,1:nvalid] <- ll_raw[i,!is.na(ll_raw[i,])]
  }
  return(ll)
}

here::i_am(".losartan_hidden_root")
here::here()
source(here::here('utils', 'losaran_r_func.R'))

renv::activate(here::here(""))
#renv::snapshot()
required_packages = c("rstan", "dplyr", "utils", "loo", "optparse", "posterior")
invisible(lapply(required_packages, require, character.only = TRUE))

run_on_cluster = 0 ## This option can be used to run the model on HPC as a `Rscript` it can automatically parse input
estimate_model = 0

if (run_on_cluster == 1) {
  # Create options to pass to script
  option_list = list(
    optparse::make_option(c('-m', '--model'), type='character', default = NULL, help = 'model', metavar = 'model'), 
    optparse::make_option(c('-f', '--flag'), type='character', default = NULL, help = 'flag', metavar = 'flag'), 
    optparse::make_option(c('-c', '--chains'), type='numeric', default = NULL, help = 'chains', metavar = 'chains'), 
    optparse::make_option(c('-i', '--iter'), type='numeric', default = NULL, help = 'iter', metavar = 'iter'), 
    optparse::make_option(c('-w', '--warmup'), type='numeric', default = NULL, help = 'warmup', metavar = 'warmp'), 
    optparse::make_option(c('-t', '--thin'), type='numeric', default = 1, help = 'thin', metavar = 'thing'), 
    optparse::make_option(c('-s', '--stimulus'), type='numeric', default = NULL, help = 'stim', metavar = 'stim'))

  # provide options in list to be callable by script
  opt_parser = optparse::OptionParser(option_list = option_list)
  opt = optparse::parse_args(opt_parser)
  model <- opt$model
  flag <- opt$flag
  chains <- opt$chains
  iter <- opt$iter
  warmup <- opt$warmup
  thin <- opt$thin
  cue <- opt$stim
  
} else {
  # Model options
  # "mm1_nh_RW_outcome"
  # "mm2_nh_RW_phase"
  # "mm3_nh_RW_outcome_phase"
  # "mm4_nh_RW_lapse"
  # "mm4_nh_RW_lapse_outcome"
  # "mm5_nh_PH_hybrid"
  # "mm5_nh_PH_hybrid_outcome"
  # "mm5_nh_PH_hybrid_kappa"
  # "mm5_nh_PH_hybrid_outcome_kappa"
  
  
  
  model <- "mm5_nh_PH_hybrid_kappa"
  # "mm5_nh_PH_hybrid_outcome_kappa"
  flag <- "nh_models" # Model family
  chains = 1 
  iter = 10
  warmup = 2
  thin = 1
  cue = 3
}


## Load prepared data 
load(here::here("data", paste0("p3_full_data_for_stanfit_cue",toString(cue),".RData")))

if (model == "mm1_nh_RW_outcome") {
  ## Model 1: RW with separate alphas for outcome
  data = list(O=O, V=V, nVis=nVis, nGroups=2, nOuts=2, gid=gid,  nSub=nSub, nObserved=nObserved, maxTr=maxTr)
  params = c("alpha", "sigma_V", "log_lik")

  } else if (model == "mm2_nh_RW_phase") {
  ## Model 1: RW with separate alphas for outcome
  data = list(O=O, V=V, nVis=nVis, nGroups=2, nPhases=2, gid=gid,  nSub=nSub, nObserved=nObserved, maxTr=maxTr)
  params = c("alpha", "sigma_V", "log_lik")
  
  } else if (model == "mm3_nh_RW_outcome_phase") {
    ## Model 1: RW with separate alphas for outcome
    data = list(O=O, V=V, nVis=nVis, nGroups=2, nPhases=2, nOuts=2, gid=gid,  nSub=nSub, nObserved=nObserved, maxTr=maxTr)
    params = c("alpha", "sigma_V", "log_lik")
    
  } else if (model == "mm4_nh_RW_lapse") {
    ## Model 1: RW with separate alphas for outcome
    data = list(O=O, V=V, nVis=nVis, nGroups=2, nOuts=2, gid=gid,  nSub=nSub, nObserved=nObserved, maxTr=maxTr)
    params = c("lapse", "alpha", "sigma_V", "log_lik")

  } else if (model == "mm5_nh_PH_hybrid_outcome_kappa") {
    ## Model 1: RW with separate alphas for outcome
    data = list(O=O, V=V, nVis=nVis, nGroups=2, nOutc=2, gid=gid,  nSub=nSub, nObserved=nObserved, maxTr=maxTr)
    params = c("eta", "alpha0", "alpha", "kappa", "sigma_V", "log_lik")
    
  } 


###################### ACTUAL FIT ##############################
modelFile = paste0(model,".stan")
print(paste0("model:", modelFile))

if (estimate_model==1) {
  print("estimating model")
  sf <- rstan::stan(file = here::here('models', modelFile),
                    data = data,pars = params,iter = iter, warmup = warmup, thin = thin, chains = chains, init = "random", algorithm = "NUTS", cores = parallel::detectCores())
  ## Save model
  sf@stanmodel@dso <- new("cxxdso")
  saveRDS(sf, file = here::here("output", paste0(modelFile, "_", flag ,"_cue",toString(cue),".rds")))
} else {
  print("loading model")
  sf = readRDS(file = here::here("output", paste0(modelFile, "_", flag ,"_cue",toString(cue),".rds")))
}


print("calculating loo and waic")
ll <- extract_log_lik(sf, merge_chains = TRUE)
loo <-loo::loo(strip_ll_of_nans(ll))
wc <-loo::waic(strip_ll_of_nans(ll))

## Also get the likelihoods separately for each group
sf <- posterior::as_draws_df(sf)

sf2 <- sf[,grepl("log_lik", colnames(sf))]
vals <- as.integer(gsub(".*,(\\d+),.*", "\\1", colnames(sf2)))
loo_gr1 <-loo::loo(strip_ll_of_nans(ll[,vals %in% which(gid==1)]))
wc_gr1 <-loo::waic(strip_ll_of_nans(ll[,vals %in% which(gid==1)]))
loo_gr2 <-loo::loo(strip_ll_of_nans(ll[,vals %in% which(gid==2)]))
wc_gr2 <-loo::waic(strip_ll_of_nans(ll[,vals %in% which(gid==2)]))

save(list=c("loo", "wc", "loo_gr1", "wc_gr1", "loo_gr2", "wc_gr2"), file= here::here("..", "output", "pl3_stuff","stan_model_fits", paste0("loo_and_waic_",modelFile,"_", flag,"_cue",toString(cue),".RData")))

print("saved loo/waic, summarizing draws")
df <- posterior::summarise_draws(sf)
write.csv(df, here::here("..", "output", "pl3_stuff", "stan_model_fits",paste0("param_summary_", modelFile, "_", flag,"_cue",toString(cue), ".csv")))
#################################################################
print("done")






