# Calibration routine using alternating years
# over the whole dataset

### load packages
library(phenor)
library(dplyr)
library(BayesianTools)

### select phenophases for calibration
pv <- "mprua65d"
classes <- "123"
burnin <- 2000
iter <- 24000

### create directory for new calibration with bayesian tools
base_dir <- file.path(
  "data/calibrations/",
  paste0("bt_run_",Sys.Date(),"_class",classes,"_",burnin,"burnin_",iter,"iter")
)

dir.create(
  base_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

### select seeds, model, years
mods <-  c("TT", "TTs", "PTT", "PTTs", "M1", "M1s")

# list driver files
driver_files <- list.files("data/pheno_net/","*35ylong.rds", full.names = TRUE)

message("Calibrating models:")
# --- run calibration for alternating years ----
message(paste0("-- calibrating: ", pv))

# read in the full driver files
# to be subset below to select the
# correct training and testing years
drivers <- readRDS(driver_files[grepl(paste0(pv,".*",classes,".*rds"), driver_files)])

# create selection criteria
# for training and testing
selection <- data.frame(
  site = drivers$site,
  year = drivers$year
) |>
  mutate(
    train = (year %% 2 == 0),
    test = (year %% 2 == 1)
  )

# select training years with stations more than length 1
selection <- selection |>
  dplyr::group_by(site) |>
  mutate(
    train = ifelse(
      length(which(train)) > 1 & train == TRUE,
      TRUE,
      FALSE
    )
  )

# select test years  with stations more than length 1
selection <- selection |>
  dplyr::group_by(site) |>
  mutate(
    test = ifelse(
      length(which(test)) > 1 & test == TRUE,
      TRUE,
      FALSE
    )
  )

# subset training and testing datasets
stn_list_train <<- pr_fm_subset(drivers, selection$train)
stn_list_test <- pr_fm_subset(drivers, selection$test)

print("fit models")

mod_out <- pr_fit_comparison(
  random_seeds = 1,
  models = mods,
  data = stn_list_train,
  method = "bayesiantools",
  control = list(
    sampler = "DEzs",
    settings = list(
      burnin = burnin,
      iterations = iter,
      nrChains = 1
    )
  )
)

geldiag <- lapply(mod_out$modelled, function(x){gelmanDiagnostics(x$opt_out, plot = F)})


### set 9999 to NA for evaluation
for(ii in 1:length(mod_out$modelled)){
  mod_out$modelled[[ii]]$predicted_values[mod_out$modelled[[ii]]$predicted_values > 1000] <- NA
}

# create output directory
dir.create(
  file.path(base_dir, pv),
  recursive = TRUE,
  showWarnings = FALSE
)

# save data as compressed RDS files
message("save as rds files")

saveRDS(
  mod_out,
  file = file.path(base_dir,pv,paste0("model_comparison_",pv,".rds")),
  compress = "xz"
)

saveRDS(
  geldiag,
  file = file.path(base_dir,pv,paste0("gelman_diagnostics_",pv,".rds")),
  compress = "xz"
)

saveRDS(
  stn_list_train,
  file = file.path(base_dir,pv,paste0("stn_list_train_",pv,".rds")),
  compress = "xz"
)

saveRDS(
  stn_list_test,
  file = file.path(base_dir,pv,paste0("stn_list_test_",pv,".rds")),
  compress = "xz"
)

