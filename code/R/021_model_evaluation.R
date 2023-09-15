##### MODEL EVALUATION #######
rm(list=ls())

### load packages
library(phenor)
library(RColorBrewer)
library(BayesianTools)

### functions
source("helpfuns.R")

### select subdir to evaluate
classes <- "123"
burnin <- 2000
iter <- 18000
date <- "2023-09-05"
subDir <- paste0("data/calibrations/bt_run_",date,"_class",classes,"_",burnin,"burnin_",iter,"iter")

### select phenophases for calibration
pv <- "mprua65d" #"mfags13d"
print(pv)

# load model parameters, testing and training data
mod_out <- readRDS(paste0(subDir,"/",pv,"/model_comparison_",pv,".rds"))
geldiag <- readRDS(paste0(subDir,"/",pv,"/gelman_diagnostics_",pv,".rds"))

stn_list_train <- readRDS(paste0(subDir,"/",pv,"/stn_list_train_",pv,".rds"))
stn_list_test <- readRDS(paste0(subDir,"/",pv,"/stn_list_test_",pv,".rds"))
train_sites <- unique(stn_list_train$site)
test_sites <- unique(stn_list_test$site)

mods <- names(mod_out$modelled)

for(sm in mods){
  pdf(file = file.path(subDir,pv,paste0("marginal_distributions_",pv,"_",sm,"_",Sys.Date(),".pdf")), width = 6, height = 6)
  par(mar = c(2,3,2,1))
  plot(mod_out$modelled[[sm]]$opt_out)
  dev.off()
}

## calculate DIC
DICvals <- sapply(mod_out$modelled, function(x) DIC(x$opt_out)$DIC)/10000
names(aicc_mat) <- mods

write.table(x = DICvals, file = file.path(subDir,pv,paste0("DIC_",pv,".txt")))

# evaluation training
eval_mods_train <- array(NA, dim=c(length(train_sites),length(mods),4))
dimnames(eval_mods_train) <- list(train_sites,mods, c("rmse","bias","pcorr","nse"))

for (mod in 1:length(mod_out$modelled)){
  for(stn in 1:length(train_sites)){
    stnind <- stn_list_train$site == train_sites[stn]
    obs <- stn_list_train$transition_dates[stnind]
    pred <- mod_out$modelled[[mod]]$predicted_values[stnind]
    eval_mods_train[stn,mod,] <- round(eval_recon(rec = pred, obs)[c("rmse","bias","pcorr","nse"),],3)
  }
}

saveRDS(
  eval_mods_train,
  file = file.path(subDir,pv,paste0("model_evaluation_train_",pv,".rds")),
  compress = "xz"
)
# evaluation testing
eval_mods_test <- array(NA, dim=c(length(test_sites),length(mods),4))
dimnames(eval_mods_test) <- list(test_sites,mods, c("rmse","bias","pcorr","nse"))

for (mod in 1:length(mod_out$modelled)){
  pred <- pr_predict(par = mod_out$modelled[[mod]]$parameters, data = stn_list_test, model = names(mod_out$modelled)[mod])
  pred[pred > 1000] <- NA ## set to NA
  
  for(stn in 1:length(test_sites)){
    stnind <- stn_list_test$site == test_sites[stn]
    obs <- stn_list_test$transition_dates[stnind]
    eval_mods_test[stn,mod,] <- round(eval_recon(rec = pred[stnind], obs)[c("rmse","bias","pcorr","nse"),],3)
  }
}

saveRDS(
  eval_mods_test,
  file = file.path(subDir,pv,paste0("model_evaluation_test_",pv,".rds")),
  compress = "xz"
)

# plot the evaluation ###
print("plot eval")

nrtest <- nrow(eval_mods_test)
nrtrain <- nrow(eval_mods_train)
nc <- length(mods)
cols <- brewer.pal(n = length(mods), "YlGnBu")

train <- eval_mods_train
test <- eval_mods_test

lwd = 0.5
png(paste0(subDir,"/",pv,"/testtrain_boxplots_",pv,"_",Sys.Date(),".png"), width = 2100, height = 900, res = 300, pointsize = 8)
par(mfcol = c(2,4), mar = c(1,2,0,1), oma = c(2,2,2,2))
rr <- c(floor(min(c(test[,,1], train[,,1]), na.rm = T)),ceiling(max(c(test[,,1], train[,,1]), na.rm = T)))
if (rr[2] > 25) rr[2] <- 25
df <- data.frame(matrix(train[,,1],nrow = nrtrain, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
mtext(side = 3,line = 0.1, "RMSE")
abline( h = seq(0,25,5), lty = 3, col = "gray75", lwd = lwd)
df <- data.frame(matrix(test[,,1],nrow = nrtest, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
axis(side = 1, at= 1:length(mods), mods)
abline( h = seq(0,25,5), lty = 3, col = "gray75", lwd = lwd)

rr <- c(floor(min(c(test[,,2], train[,,2]), na.rm = T)),ceiling(max(c(test[,,2], train[,,2]), na.rm = T)))
df <- data.frame(matrix(train[,,2],nrow = nrtrain, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
mtext(side = 3,line = 0.1, "Mean bias")
abline( h = seq(-25,25,5), lty = 3, col = "gray75", lwd = lwd)
df <- data.frame(matrix(test[,,2],nrow = nrtest, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
axis(side = 1, at = 1:length(mods), mods)
abline( h = seq(-25,25,5), lty = 3, col = "gray75", lwd = lwd)

rr <- c(floor(min(c(test[,,3], train[,,3]), na.rm = T)),ceiling(max(c(test[,,3], train[,,3]), na.rm = T)))
df <- data.frame(matrix(train[,,3],nrow = nrtrain, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
mtext(side = 3,line = 0.1, "Pearson correlation")
abline( h = seq(-1,1,0.2), lty = 3, col = "gray75", lwd = lwd)
df <- data.frame(matrix(test[,,3],nrow = nrtest, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
axis(side = 1, at = 1:length(mods), mods)
abline( h = seq(-1,1,0.2), lty = 3, col = "gray75", lwd = lwd)

rr <- c(floor(min(c(test[,,4], train[,,4]), na.rm = T)),ceiling(max(c(test[,,4], train[,,4]), na.rm = T)))
if (rr[1] < -3) rr[1] <- -3
df <- data.frame(matrix(train[,,4],nrow = nrtrain, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
mtext(side = 3,line = 0.1, "NSE")
abline( h = seq(-4,1,0.5), lty = 3, col = "gray75", lwd = lwd)
df <- data.frame(matrix(test[,,4],nrow = nrtest, ncol = nc, byrow = F))
boxplot(df, ylim = rr, xaxt = "n", col = cols)
axis(side = 1, at = 1:length(mods), mods)
abline( h = seq(-8,1,0.5), lty = 3, col = "gray75", lwd = lwd)
dev.off()


